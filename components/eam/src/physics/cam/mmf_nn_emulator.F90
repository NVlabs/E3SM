module mmf_nn_emulator

use constituents,    only: pcnst
use shr_kind_mod,    only: r8 => shr_kind_r8
use ppgrid,          only: pcols, pver, pverp
use cam_history,     only: outfld, addfld, horiz_only
use physconst,       only: gravit,cpair,latvap,latice
use spmd_utils,      only: masterproc,iam
use camsrfexch,      only: cam_out_t, cam_in_t
use constituents,    only: cnst_get_ind
use physics_types,   only: physics_state,  physics_ptend, physics_ptend_init
use cam_logfile,     only: iulog
use physics_buffer,  only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
use cam_history_support, only: pflds, fieldname_lenp2
use cam_abortutils,  only: endrun
use string_utils,    only: to_lower
use phys_grid,       only: get_lat_p, get_lon_p, get_rlat_p, get_rlon_p

!-------- torch fortran --------

use torch_ftn
use iso_fortran_env
!--------------------------------------

  implicit none
  save 

  private
  ! Define variables for this entire module
  ! These are nameliest variables.
  ! If not specified in atm_in, the following defaults values are used.
  integer :: inputlength  = 1525     ! length of NN input vector
  integer :: outputlength = 368     ! length of NN output vector
  logical :: input_rh     = .false.  ! toggle to switch from q --> RH input
  integer :: cb_spinup_step = 72 
  logical :: cb_use_input_prectm1 = .false.  ! use previous timestep PRECT for input variable 
  logical :: cb_partial_coupling  = .false.
  character(len=fieldname_lenp2) :: cb_partial_coupling_vars(pflds)
  character(len=256) :: cb_nn_var_combo = 'v4' ! nickname for a specific NN in/out variable combination, currently support v2 or v4
  character(len=256)    :: cb_torch_model   ! absolute filepath for a torchscript model txt file

  ! option to mix the NN output with an SP prediction with customized partition scheduling
  !tendency sent to GCM grid would be: ratio * NN output + (1-ratio) * SP output
  logical :: cb_do_ramp = .false. ! option to turn on customized ratio scheduling
  real(r8) :: cb_ramp_factor = 1.0
  character(len=256)    :: cb_ramp_option = 'constant' ! 'constant', 'linear', 'step'
  ! if cb_ramp_option = 'constant', ratio = cb_ramp_factor constantly in time
  ! if cb_ramp_option = 'linear', do linearly ramping up the ratio value
  !     ratio = cb_ramp_factor * (time/cb_ramp_linear_steps) when time < cb_ramp_linear_steps, ratio = cb_ramp_factor when time >= cb_ramp_linear_steps
  ! if cb_ramp_option = 'step', periodically switch the ratio value between cb_ramp_factor and 0.0
  !     ratio = cb_ramp_factor for cb_ramp_step_1steps, then 0.0 for cb_ramp_step_0steps, then back to cb_ramp_factor
  integer :: cb_ramp_linear_steps = 0
  integer :: cb_ramp_step_0steps = 10
  integer :: cb_ramp_step_1steps = 20
  
  logical :: cb_strato_water_constraint = .false. ! zero out cloud (qc and qi) above tropopause in the NN output
  real(r8) :: dtheta_thred = 10.0 ! threshold for determine the tropopause (currently is p<400hPa and dtheta/dz>dtheta_thred = 10K/km) 

  type(torch_module), allocatable :: torch_mod(:)

  ! local
  logical :: cb_top_levels_zero_out = .true.
  integer :: cb_n_levels_zero = 12 ! top n levels to zero out

#ifdef MMF_NN_EMULATOR
  public neural_net, init_neural_net, mmf_nn_emulator_readnl, &
         cb_partial_coupling, cb_partial_coupling_vars, cb_spinup_step, cb_do_ramp, cb_ramp_linear_steps, cb_ramp_option, cb_ramp_factor, cb_ramp_step_0steps, cb_ramp_step_1steps
#else
  public mmf_nn_emulator_readnl, cb_partial_coupling, cb_partial_coupling_vars, cb_spinup_step, cb_do_ramp, cb_ramp_linear_steps, cb_ramp_option, cb_ramp_factor, cb_ramp_step_0steps, cb_ramp_step_1steps
#endif
  
contains

#ifdef MMF_NN_EMULATOR
  subroutine neural_net (ptend, state, state_aphys1, pbuf, cam_in, cam_out, coszrs, solin, ztodt, lchnk)
 ! note state is meant to have the "BP" state saved earlier. 

   implicit none

   type(physics_ptend),intent(out)    :: ptend 
   type(physics_state), intent(in)    :: state
   type(physics_state), intent(in)    :: state_aphys1
   type(physics_buffer_desc), pointer :: pbuf(:)  
   type(cam_in_t),intent(in)          :: cam_in
   type(cam_out_t),     intent(inout) :: cam_out 
   real(r8), intent(in)               :: coszrs(pcols)
   real(r8), intent(in)               :: solin(pcols)
   real(r8), intent(in)               :: ztodt
   integer, intent(in)                :: lchnk


   ! local variables
   real(r8) :: input(pcols,inputlength)
   real(r8) :: output(pcols,outputlength)
   integer :: i,k,ncol,ixcldice,ixcldliq,ii,kk,idx_trop(pcols),kens, j
   real (r8) :: s_bctend(pcols,pver), q_bctend(pcols,pver), qc_bctend(pcols,pver), qi_bctend(pcols,pver), qafter, safter, &
                u_bctend(pcols,pver), v_bctend(pcols,pver)
   logical :: do_constraints
   logical :: lq(pcnst)
   real(r8) ::rh_loc
   integer :: prec_dp_idx, snow_dp_idx
   ! convective precipitation variables
   real(r8), pointer :: prec_dp(:)                ! total precip rate [m/s]
   real(r8), pointer :: snow_dp(:)                ! snow precip rate [m/s]

   real(r8), pointer, dimension(:)   :: lhflx, shflx, taux, tauy ! (/pcols/)
   real(r8), pointer, dimension(:,:) :: ozone, ch4, n2o ! (/pcols,pver/)

  !  type(torch_module) :: torch_mod
   type(torch_tensor_wrap) :: input_tensors
   type(torch_tensor) :: out_tensor
   real(real32) :: input_torch(inputlength, pcols)
   real(real32), pointer :: output_torch(:, :)
   real(r8) :: math_pi

   math_pi = 3.14159265358979323846_r8

   ncol  = state%ncol
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   lq(:)        = .FALSE.
   lq(1)        = .TRUE. ! water vapor
   lq(ixcldliq) = .TRUE. ! cloud liquid
   lq(ixcldice) = .TRUE. ! cloud ice
   call physics_ptend_init(ptend, state%psetcols, 'neural-net', & ! Initialize local physics_ptend object
                           ls=.true., &
                           lq=lq,     &
                           lu=.true., &
                           lv=.true.  &
                          )

   do_constraints = .true.
   
   s_bctend(:,:)  = 0.
   q_bctend(:,:)  = 0.
   qc_bctend(:,:) = 0.
   qi_bctend(:,:) = 0.
   u_bctend(:,:)  = 0.
   v_bctend(:,:)  = 0.

   ! Associate pointers with pbuf fields
   call pbuf_get_field(pbuf, pbuf_get_index('LHFLX'), lhflx)
   call pbuf_get_field(pbuf, pbuf_get_index('SHFLX'), shflx)
   call pbuf_get_field(pbuf, pbuf_get_index('TAUX'),  taux)
   call pbuf_get_field(pbuf, pbuf_get_index('TAUY'),  tauy)
   call pbuf_get_field(pbuf, pbuf_get_index('ozone'), ozone)
   call pbuf_get_field(pbuf, pbuf_get_index('CH4'),   ch4)
   call pbuf_get_field(pbuf, pbuf_get_index('N2O'),   n2o)

! v2 input variables:
! ['state_t', 'state_q0001', 'state_q0002', 'state_q0003', 'state_u', 'state_v', 'state_ps', 'pbuf_SOLIN', 'pbuf_LHFLX', 'pbuf_SHFLX', 'pbuf_TAUX', 'pbuf_TAUY', 'pbuf_COSZRS', 'cam_in_ALDIF', 'cam_in_ALDIR', 'cam_in_ASDIF', 'cam_in_ASDIR', 'cam_in_LWUP', 'cam_in_ICEFRAC', 'cam_in_LANDFRAC', 'cam_in_OCNFRAC', 'cam_in_SNOWHICE', 'cam_in_SNOWHLAND', 'pbuf_ozone', 'pbuf_CH4', 'pbuf_N2O']
! v2 output variables:
! ['ptend_t', 'ptend_q0001', 'ptend_q0002', 'ptend_q0003', 'ptend_u', 'ptend_v', 'cam_out_NETSW', 'cam_out_FLWDS', 'cam_out_PRECSC', 'cam_out_PRECC', 'cam_out_SOLS', 'cam_out_SOLL', 'cam_out_SOLSD', 'cam_out_SOLLD']

select case (to_lower(trim(cb_nn_var_combo)))

    case('v2')
      input(:ncol,0*pver+1:1*pver) = state%t(1:ncol,1:pver)          ! state_t
      input(:ncol,1*pver+1:2*pver) = state%q(1:ncol,1:pver,1)        ! state_q0001
      input(:ncol,2*pver+1:3*pver) = state%q(1:ncol,1:pver,ixcldliq)
      input(:ncol,3*pver+1:4*pver) = state%q(1:ncol,1:pver,ixcldice)
      input(:ncol,4*pver+1:5*pver) = state%u(1:ncol,1:pver)          ! state_u
      input(:ncol,5*pver+1:6*pver) = state%v(1:ncol,1:pver)          ! state_v
      !gas
      input(:ncol,6*pver+1:7*pver) = ozone(:ncol,1:pver)            ! pbuf_ozone
      input(:ncol,7*pver+1:8*pver) = ch4(:ncol,1:pver)             ! pbuf_CH4
      input(:ncol,8*pver+1:9*pver) = n2o(:ncol,1:pver)             ! pbuf_N2O
      ! 2d vars e.g., ps, solin
      input(:ncol,9*pver+1) = state%ps(1:ncol)                      ! state_ps
      input(:ncol,9*pver+2) = solin(1:ncol)                         ! pbuf_SOLIN
      input(:ncol,9*pver+3) = lhflx(1:ncol)                         ! pbuf_LHFLX
      input(:ncol,9*pver+4) = shflx(1:ncol)                         ! pbuf_SHFLX
      input(:ncol,9*pver+5) = taux(1:ncol)                          ! pbuf_TAUX
      input(:ncol,9*pver+6) = tauy(1:ncol)                          ! pbuf_TAUY
      input(:ncol,9*pver+7) = coszrs(1:ncol)                        ! pbuf_COSZRS
      input(:ncol,9*pver+8) = cam_in%ALDIF(:ncol)                   ! cam_in_ALDIF
      input(:ncol,9*pver+9) = cam_in%ALDIR(:ncol)                   ! cam_in_ALDIR
      input(:ncol,9*pver+10) = cam_in%ASDIF(:ncol)                  ! cam_in_ASDIF
      input(:ncol,9*pver+11) = cam_in%ASDIR(:ncol)                  ! cam_in_ASDIR
      input(:ncol,9*pver+12) = cam_in%LWUP(:ncol)                   ! cam_in_LWUP
      input(:ncol,9*pver+13) = cam_in%ICEFRAC(:ncol)                ! cam_in_ICEFRAC
      input(:ncol,9*pver+14) = cam_in%LANDFRAC(:ncol)               ! cam_in_LANDFRAC
      input(:ncol,9*pver+15) = cam_in%OCNFRAC(:ncol)                ! cam_in_OCNFRAC
      input(:ncol,9*pver+16) = cam_in%SNOWHICE(:ncol)               ! cam_in_SNOWHICE
      input(:ncol,9*pver+17) = cam_in%SNOWHLAND(:ncol)              ! cam_in_SNOWHLAND
      ! RH conversion
      if (input_rh) then ! relative humidity conversion for input
        do i = 1,ncol
          do k=1,pver
            ! Port of tom's RH =  Rv*p*qv/(R*esat(T))
            rh_loc = 461.*state%pmid(i,k)*state%q(i,k,1)/(287.*tom_esat(state%t(i,k))) ! note function tom_esat below refercing SAM's sat.F90
#ifdef RHDEBUG
            if (masterproc) then
              write (iulog,*) 'RHDEBUG:p,q,T,RH=',state%pmid(i,k),state%q(i,k,1),state%t(i,k),rh_loc
            endif
#endif
            input(i,1*pver+k) = rh_loc
          end do
        end do
      end if

    case('v4')
      input(:ncol,0*pver+1:1*pver) = state%t(1:ncol,1:pver)          ! state_t
      input(:ncol,1*pver+1:2*pver) = state%q(1:ncol,1:pver,1)        ! state_q0001
      input(:ncol,2*pver+1:3*pver) = state%q(1:ncol,1:pver,ixcldliq)
      input(:ncol,3*pver+1:4*pver) = state%q(1:ncol,1:pver,ixcldice)
      input(:ncol,4*pver+1:5*pver) = state%u(1:ncol,1:pver)          ! state_u
      input(:ncol,5*pver+1:6*pver) = state%v(1:ncol,1:pver)          ! state_v
      ! state dynamics tendencies (large-scale forcings)
      input(:ncol,6*pver+1:7*pver) = (state%t(1:ncol,1:pver)-state_aphys1%t(1:ncol,1:pver))/1200. ! state_t_dyn
      ! state_q0_dyn
      input(:ncol,7*pver+1:8*pver) = (state%q(1:ncol,1:pver,1)-state_aphys1%q(1:ncol,1:pver,1) + state%q(1:ncol,1:pver,ixcldliq)-state_aphys1%q(1:ncol,1:pver,ixcldliq) + state%q(1:ncol,1:pver,ixcldice)-state_aphys1%q(1:ncol,1:pver,ixcldice))/1200.
      input(:ncol,8*pver+1:9*pver) = (state%u(1:ncol,1:pver)-state_aphys1%u(1:ncol,1:pver))/1200. ! state_u_dyn
      input(:ncol,9*pver+1:10*pver) = state%t_adv(2,1:ncol,1:pver)
      input(:ncol,10*pver+1:11*pver) = state%q_adv(2,1:ncol,1:pver,1) + state%q_adv(2,1:ncol,1:pver,ixcldliq) + state%q_adv(2,1:ncol,1:pver,ixcldice)
      input(:ncol,11*pver+1:12*pver) = state%u_adv(2,1:ncol,1:pver)
      ! previous state physics tendencies
      input(:ncol,12*pver+1:13*pver) = state%t_phy(1,1:ncol,1:pver)
      input(:ncol,13*pver+1:14*pver) = state%q_phy(1,1:ncol,1:pver,1)
      input(:ncol,14*pver+1:15*pver) = state%q_phy(1,1:ncol,1:pver,ixcldliq)
      input(:ncol,15*pver+1:16*pver) = state%q_phy(1,1:ncol,1:pver,ixcldice)
      input(:ncol,16*pver+1:17*pver) = state%u_phy(1,1:ncol,1:pver)
      ! 2-step in the past physics tendencies
      input(:ncol,17*pver+1:18*pver) = state%t_phy(2,1:ncol,1:pver)
      input(:ncol,18*pver+1:19*pver) = state%q_phy(2,1:ncol,1:pver,1)
      input(:ncol,19*pver+1:20*pver) = state%q_phy(2,1:ncol,1:pver,ixcldliq)
      input(:ncol,20*pver+1:21*pver) = state%q_phy(2,1:ncol,1:pver,ixcldice)
      input(:ncol,21*pver+1:22*pver) = state%u_phy(2,1:ncol,1:pver)
      !gas
      input(:ncol,22*pver+1:23*pver) = ozone(:ncol,1:pver)            ! pbuf_ozone
      input(:ncol,23*pver+1:24*pver) = ch4(:ncol,1:pver)             ! pbuf_CH4
      input(:ncol,24*pver+1:25*pver) = n2o(:ncol,1:pver)             ! pbuf_N2O
      ! 2d vars e.g., ps, solin
      input(:ncol,25*pver+1) = state%ps(1:ncol)                      ! state_ps
      input(:ncol,25*pver+2) = solin(1:ncol)                         ! pbuf_SOLIN
      input(:ncol,25*pver+3) = lhflx(1:ncol)                         ! pbuf_LHFLX
      input(:ncol,25*pver+4) = shflx(1:ncol)                         ! pbuf_SHFLX
      input(:ncol,25*pver+5) = taux(1:ncol)                          ! pbuf_TAUX
      input(:ncol,25*pver+6) = tauy(1:ncol)                          ! pbuf_TAUY
      input(:ncol,25*pver+7) = coszrs(1:ncol)                        ! pbuf_COSZRS
      input(:ncol,25*pver+8) = cam_in%ALDIF(:ncol)                   ! cam_in_ALDIF
      input(:ncol,25*pver+9) = cam_in%ALDIR(:ncol)                   ! cam_in_ALDIR
      input(:ncol,25*pver+10) = cam_in%ASDIF(:ncol)                  ! cam_in_ASDIF
      input(:ncol,25*pver+11) = cam_in%ASDIR(:ncol)                  ! cam_in_ASDIR
      input(:ncol,25*pver+12) = cam_in%LWUP(:ncol)                   ! cam_in_LWUP
      input(:ncol,25*pver+13) = cam_in%ICEFRAC(:ncol)                ! cam_in_ICEFRAC
      input(:ncol,25*pver+14) = cam_in%LANDFRAC(:ncol)               ! cam_in_LANDFRAC
      input(:ncol,25*pver+15) = cam_in%OCNFRAC(:ncol)                ! cam_in_OCNFRAC
      input(:ncol,25*pver+16) = cam_in%SNOWHICE(:ncol)               ! cam_in_SNOWHICE
      input(:ncol,25*pver+17) = cam_in%SNOWHLAND(:ncol)              ! cam_in_SNOWHLAND
      !5 placeholder for future input (which we plan to use some 2d vars at t-1, but in the end did not implement)
      input(:ncol,25*pver+18:25*pver+22) = 0._r8
      ! cos lat and sin lat
      do i = 1,ncol ! lat is get_lat_p(lchnk,i), 23/24 needs cos/sin
        input(i,25*pver+23) = cos(get_rlat_p(lchnk,i))
        input(i,25*pver+24) = sin(get_rlat_p(lchnk,i))
      end do
      input(:ncol,25*pver+25) = 0._r8               ! icol ! can be 1-384 in future but not implemented yet
      ! RH conversion
      if (input_rh) then ! relative humidity conversion for input
         do i = 1,ncol
           do k=1,pver
             ! Port of tom's RH =  Rv*p*qv/(R*esat(T))
             rh_loc = 461.*state%pmid(i,k)*state%q(i,k,1)/(287.*tom_esat(state%t(i,k))) ! note function tom_esat below refercing SAM's sat.F90
             input(i,1*pver+k) = rh_loc
           end do
         end do
      end if
end select


    ! do the torch inference    
    input_torch(:,:) = 0.
    do i=1,ncol
      do k=1,inputlength
        input_torch(k,i) = input(i,k)
      end do
    end do
    !print *, "Creating input tensor"
    call input_tensors%create
    !print *, "Adding input data"
    call input_tensors%add_array(input_torch)
    call torch_mod(1)%forward(input_tensors, out_tensor, flags=module_use_inference_mode)
    call out_tensor%to_array(output_torch)
    do i=1, ncol
      do k=1,outputlength
        output(i,k) = output_torch(k,i)
      end do
    end do


    !!! do post processing some non-negative constraints + option to do stratosphere water denial

    ! Manually applying ReLU activation for positive-definite variables
    ! [TODO] for ensemble, ReLU should be moved before ens-averaging
    do i=1,ncol
      ! ReLU for the last 8 variables 
      do k=outputlength-7,outputlength
        output(i,k) = max(output(i,k), 0.)
      end do
      ! ! tiny flwds
      ! k=6*pver+2
      ! output(i,k) = max(output(i,k), tiny(output(i,k))) ! flwds
      !                                                   ! preventing flwds==0 error

      ! zero out surface solar fluxes when local time is at night
      if (coszrs(i) .le. 0.) then
        output(i,6*pver+1) = 0. ! netsw
        output(i,6*pver+5) = 0. ! sols
        output(i,6*pver+6) = 0. ! soll
        output(i,6*pver+7) = 0. ! solsd
        output(i,6*pver+8) = 0. ! solld
      endif
    end do
 
 
 ! ---------- 1. NN output to atmosphere forcing --------
 ! ['TBCTEND', 'QBCTEND','CLDLIQBCTEND','CLDICEBCTEND']
    s_bctend (:ncol,1:pver) = output(1:ncol,0*pver+1:1*pver)*cpair ! K/s --> J/kg/s (ptend expects that)
    q_bctend (:ncol,1:pver) = output(1:ncol,1*pver+1:2*pver)       ! kg/kg/s
    qc_bctend(:ncol,1:pver) = output(1:ncol,2*pver+1:3*pver)       ! kg/kg/s
    qi_bctend(:ncol,1:pver) = output(1:ncol,3*pver+1:4*pver)       ! kg/kg/s
    u_bctend (:ncol,1:pver) = output(1:ncol,4*pver+1:5*pver)       ! m/s/s
    v_bctend (:ncol,1:pver) = output(1:ncol,5*pver+1:6*pver)       ! m/s/s
 
 ! deny any moisture activity and remove all clouds in the stratosphere:
    if (cb_strato_water_constraint) then
     do i=1,ncol
       call detect_tropopause(state%t(i,:),state%exner(i,:),state%zm(i,:),state%pmid(i,:),idx_trop(i))
       q_bctend (i,1:idx_trop(i)) = 0.
       qc_bctend(i,1:idx_trop(i)) = -state%q(i,1:idx_trop(i),ixcldliq)/ztodt
       qi_bctend(i,1:idx_trop(i)) = -state%q(i,1:idx_trop(i),ixcldice)/ztodt
     end do
     call outfld('TROP_IND', idx_trop(:ncol)*1._r8, ncol, state%lchnk)
    end if
 
 ! -- atmos positivity constraints ---- 
    if (do_constraints) then
    do i=1,ncol
      do k=1,pver
 
 ! energy positivity:
        safter = state%s(i,k) + s_bctend(i,k)*ztodt ! predicted DSE after NN tendency
        if (safter .lt. 0.) then ! can only happen when bctend < 0...
          s_bctend(i,k) = s_bctend(i,k) + abs(safter)/ztodt ! in which case reduce cooling rate
          write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute temperature, corrected but BEWARE!!!'
          write (iulog,*) '' ! [TODO] printout lat/lon and error magnitude
        endif
 
  ! vapor positivity:
        qafter = state%q(i,k,1) + q_bctend(i,k)*ztodt ! predicted vapor after NN tendency
        if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
          q_bctend(i,k) = q_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
          write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute q, corrected but BEWARE!!!'
        endif
 
  ! liquid positivity:
        qafter = state%q(i,k,ixcldliq) + qc_bctend(i,k)*ztodt ! predicted liquid after NN tendency
        if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
          qc_bctend(i,k) = qc_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
          !write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute qc, corrected but BEWARE!!!'
        endif
 ! ice positivity:
        qafter = state%q(i,k,ixcldice) + qi_bctend(i,k)*ztodt ! predicted ice after NN tendency
        if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
          qi_bctend(i,k) = qi_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
          !write (iulog,*) 'HEY MMF_NN_EMULATOR made a negative absolute qi, corrected but BEWARE!!!'
        endif
      end do
    end do
    endif
 ! Wire to ptend:
     ptend%s(:ncol,:pver)          = s_bctend(:ncol,:pver)
     ptend%q(:ncol,:pver,1)        = q_bctend(:ncol,:pver)
     ptend%q(:ncol,:pver,ixcldliq) = qc_bctend(:ncol,:pver)
     ptend%q(:ncol,:pver,ixcldice) = qi_bctend(:ncol,:pver)
     ptend%u(:ncol,:pver)          = u_bctend(:ncol,:pver)
     ptend%v(:ncol,:pver)          = v_bctend(:ncol,:pver)
 
 ! for q{1,2,3}, u, v tendencies, top levels [1 to 12] are not coupled
 ! no std dev except ptend%q1[1:3], whose magnitudes are smaller that other levels by orders of magnitude
     if (cb_top_levels_zero_out) then
       ptend%q(:ncol,1:cb_n_levels_zero,1)        = 0. 
       ptend%q(:ncol,1:cb_n_levels_zero,ixcldliq) = 0. 
       ptend%q(:ncol,1:cb_n_levels_zero,ixcldice) = 0. 
       ptend%u(:ncol,1:cb_n_levels_zero)          = 0. 
       ptend%v(:ncol,1:cb_n_levels_zero)          = 0. 
     end if
 
 ! ------------- 2. NN output to land forcing ---------
 !!! Sungduk: It works, but I wrote it again to add 'ocean only coupling' option
 !!!          (#MMF_NN_EMULATOR_OCN_ONLY)
 !!! 
 !!!    ! These are the cam_out members that are assigned in cam_export: prec_dp, snow_dp,
 !!!    ! and so saved to pbuf, instead.
 !!!    ! SY: Note that this uses SPCAM's or CAM's pbuf register: crm_physics_register(), convect_deep_register()
 !!!    !     Once, SP is entirely removed, cam_out%prec{c,l,sc,sl} should be directly updated here.
 
    prec_dp_idx = pbuf_get_index('PREC_DP', errcode=i) ! Query physics buffer index
    snow_dp_idx = pbuf_get_index('SNOW_DP', errcode=i)
    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp) ! Associate pointers withphysics buffer fields
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
 
    do i = 1,ncol
 ! SY: debugging
 !     allowing surface coupling over ocean only
#ifdef MMF_NN_EMULATOR_OCN_ONLY 
      if (cam_in%ocnfrac(i) .eq. 1.0_r8) then
#endif
        cam_out%netsw(i) = output(i,6*pver+1)
        cam_out%flwds(i) = output(i,6*pver+2)
        snow_dp(i)       = output(i,6*pver+3)
        prec_dp(i)       = output(i,6*pver+4)
        cam_out%sols(i)  = output(i,6*pver+5)
        cam_out%soll(i)  = output(i,6*pver+6)
        cam_out%solsd(i) = output(i,6*pver+7)
        cam_out%solld(i) = output(i,6*pver+8)
#ifdef MMF_NN_EMULATOR_OCN_ONLY
      end if
#endif
    end do 


end subroutine neural_net
#endif

#ifdef MMF_NN_EMULATOR
  subroutine init_neural_net()

    implicit none

    integer :: i, k

    allocate(torch_mod (1))
    call torch_mod(1)%load(trim(cb_torch_model), 0) !0 is not using gpu, for now just use cpu for NN inference
    !call torch_mod(1)%load(trim(cb_torch_model), module_use_device) will use gpu if available
    
  ! add diagnostic output fileds
  call addfld ('TROP_IND',horiz_only,   'A', '1', 'lev index for tropopause')

  end subroutine init_neural_net
#endif

  real(r8) function tom_esat(T) 
  ! For consistency with the python port of Tom's RH-calculator, this is how it
  ! was done in the training environment (Caution: could be porting bugs here)
    implicit none
    real(r8) T
    real(r8), parameter :: T0 = 273.16
    real(r8), parameter :: T00 = 253.16
    real(r8), external :: esatw_crm,esati_crm ! register functions from crm source.
    real(r8) :: omtmp,omega
    omtmp = (T-T00)/(T0-T00)
    omega = max(0.,min(1.,omtmp))
!tf.where(T>T0,eliq(T),tf.where(T<T00,eice(T),(omega*eliq(T)+(1-omega)*eice(T))))
    if (T .gt. T0) then
      tom_esat = tom_eliq(T)
    elseif (T .lt. T00) then
      tom_esat = tom_eice(T)
    else
      tom_esat = omega*tom_eliq(T) + (1.-omega)*tom_eice(T)
    endif
  end

  real(r8) function tom_eliq(T)
    implicit none
    real(r8) T
    real(r8), parameter :: T0 = 273.16
    real(r8), parameter :: cliq = -80. 
    real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
       6.11239921, 0.443987641, 0.142986287e-1, &
       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real(r8) :: dt
    dt = max(cliq,T-T0)
    tom_eliq = 100.*(a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))  
  end 


  real(r8) function tom_eice(T)
    implicit none
    real(r8) T
    real(r8), parameter :: T0 = 273.16
    real(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/       
    real(r8) cice(6)
    real(r8) dt
    dt = T-T0
    cice(1) = 273.15
    cice(2) = 185.
    cice(3) = -100.
    cice(4) = 0.00763685
    cice(5) = 0.000151069
    cice(6) = 7.48215e-07
    if (T .gt. cice(1)) then
      tom_eice = tom_eliq(T)
    else if (T .le. cice(2)) then
      tom_eice = 100.*(cice(4) + max(cice(3),dt)*(cice(5)+max(cice(3),dt)*cice(6))) 
    else
      tom_eice = 100.*(a0 +dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))))
    end if
  end      


 subroutine detect_tropopause (t,exner,zmid,pmid,klev_crmtop)
   real(r8), intent(in) :: t(pver),exner(pver),zmid(pver),pmid(pver)
   integer, intent(out) :: klev_crmtop
   integer :: k
   real (r8) :: theta(pver),dthetadz

   do k=1,pver
     theta(k) = t(k)*exner(k)
   end do

   klev_crmtop = 1

   do k=2,pver-1
     dthetadz = (theta(k-1)-theta(k+1))/(zmid(k-1)-zmid(k+1))*1000. ! K/km
     ! assume theta in K and pmid in Pa then
     if (pmid(k) .le. 40000. .and. dthetadz > dtheta_thred) then
       klev_crmtop = k
     endif
   end do
  end subroutine detect_tropopause

  ! Read namelist variables.
  subroutine mmf_nn_emulator_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr, f
      character(len=*), parameter :: subname = 'mmf_nn_emulator_readnl'
      
      namelist /mmf_nn_emulator_nl/ inputlength, outputlength, input_rh, &
                           cb_partial_coupling, cb_partial_coupling_vars,&
                           cb_use_input_prectm1, &
                           cb_nn_var_combo, &
                           cb_torch_model, cb_spinup_step, &
                           cb_do_ramp, cb_ramp_linear_steps, &
                           cb_ramp_option, cb_ramp_factor, cb_ramp_step_0steps, cb_ramp_step_1steps, &
                           cb_strato_water_constraint, dtheta_thred

      ! Initialize 'cb_partial_coupling_vars'
      do f = 1, pflds
        cb_partial_coupling_vars(f) = ' '
      end do

      ierr = 0
      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'mmf_nn_emulator_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, mmf_nn_emulator_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(inputlength,  1,                 mpiint,  0, mpicom)
      call mpibcast(outputlength, 1,                 mpiint,  0, mpicom)
      call mpibcast(input_rh,     1,                 mpilog,  0, mpicom)
      call mpibcast(cb_use_input_prectm1,1,          mpilog,  0, mpicom)
      call mpibcast(cb_partial_coupling, 1,          mpilog,  0, mpicom)
      call mpibcast(cb_partial_coupling_vars, len(cb_partial_coupling_vars(1))*pflds, mpichar, 0, mpicom)
      call mpibcast(cb_nn_var_combo, len(cb_nn_var_combo), mpichar,  0, mpicom)
      call mpibcast(cb_torch_model, len(cb_torch_model), mpichar, 0, mpicom)
      call mpibcast(cb_spinup_step,   1, mpiint,  0, mpicom)
      call mpibcast(cb_do_ramp, 1,                  mpilog,  0, mpicom)
      call mpibcast(cb_ramp_linear_steps, 1,            mpiint,  0, mpicom)
      call mpibcast(cb_ramp_option, len(cb_ramp_option), mpichar,  0, mpicom)
      call mpibcast(cb_ramp_factor, 1,            mpir8,  0, mpicom)
      call mpibcast(cb_ramp_step_0steps, 1,            mpiint,  0, mpicom)
      call mpibcast(cb_ramp_step_1steps, 1,            mpiint,  0, mpicom)
      call mpibcast(cb_strato_water_constraint,   1, mpilog,  0, mpicom)
      call mpibcast(dtheta_thred, 1,            mpir8,  0, mpicom)


      ! [TODO] check ierr for each mpibcast call
      ! if (ierr /= 0) then
      !    call endrun(subname // ':: ERROR broadcasting namelist variable cb_partial_coupling_vars')
      ! end if
#endif

   end subroutine mmf_nn_emulator_readnl

  function shuffle_1d(array_1d) result(array_shuffled)
  ! Shuffling the entries of 1-d INTEGER array
  ! (using the Knuth shuffle algorithm: https://en.wikipedia.org/wiki/Fisher–Yates_shuffle)
    implicit none
    integer, intent(in) :: array_1d(:)
    integer :: array_shuffled(size(array_1d))
    integer :: j, k, tmp
    real :: u

    array_shuffled(:) = array_1d(:)

    j    = size(array_1d)
    do while (j > 1)
       call random_seed
       call random_number(u)
       k = 1 + FLOOR(j * u)
       tmp = array_shuffled(j)
       array_shuffled(j) = array_shuffled(k)
       array_shuffled(k) = tmp
       j = j -1
    end do
  end function shuffle_1d

end module mmf_nn_emulator
