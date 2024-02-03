module climsim

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

!-------- NEURAL-FORTRAN (FKB) --------
! use mod_kinds, only: ik, rk
use mod_network, only: network_type
!--------------------------------------

  implicit none
  save 

  private
  ! Define variables for this entire module
  ! These are nameliest variables.
  ! If not specified in atm_in, the following defaults values are used.
  integer :: inputlength  = 425     ! length of NN input vector
  integer :: outputlength = 368     ! length of NN output vector
  logical :: input_rh     = .false.  ! toggle to switch from q --> RH input
  logical :: cb_use_input_prectm1 = .false.  ! use previous timestep PRECT for input variable 
  character(len=256)    :: cb_fkb_model   ! absolute filepath for a fkb model txt file
  character(len=256)    :: cb_inp_sub     ! absolute filepath for a inp_sub.txt
  character(len=256)    :: cb_inp_div     ! absolute filepath for a inp_div.txt
  character(len=256)    :: cb_out_scale   ! absolute filepath for a out_scale.txt
  logical :: cb_partial_coupling  = .false.
  character(len=fieldname_lenp2) :: cb_partial_coupling_vars(pflds)
  character(len=256) :: cb_nn_var_combo = 'v2' ! nickname for a specific NN in/out variable combination

  type(network_type), allocatable :: climsim_net(:)
  real(r8), allocatable :: inp_sub(:)
  real(r8), allocatable :: inp_div(:)
  real(r8), allocatable :: out_scale(:)

  logical :: cb_do_ensemble  = .false. ! ensemble model inference
  integer :: cb_ens_size               ! # of ensemble models
  integer :: max_nn_ens = 100 ! Max. ensemble size is arbitrarily set to 100.
  character(len=256), allocatable :: cb_ens_fkb_model_list(:) ! absolute filepath for models

  integer :: cb_random_ens_size = 0 ! For stochastic ensemble,
                                    ! # of ensemble subset to be averaged

  ! local
  integer, allocatable :: ens_ind_shuffled(:)
  logical :: cb_do_random_ensemble = .false.
  logical :: cb_top_levels_zero_out = .true.
  integer :: cb_n_levels_zero = 12 ! top n levels to zero out

  public neural_net, init_neural_net, climsim_readnl, &
         cb_partial_coupling, cb_partial_coupling_vars
  
contains

  subroutine neural_net (ptend, state, pbuf, cam_in, cam_out, coszrs, solin, ztodt)
 ! note state is meant to have the "BP" state saved earlier. 

   implicit none

   type(physics_ptend),intent(out)    :: ptend 
   type(physics_state), intent(in)    :: state
   type(physics_buffer_desc), pointer :: pbuf(:)  
   type(cam_in_t),intent(in)          :: cam_in
   type(cam_out_t),     intent(inout) :: cam_out 
   real(r8), intent(in)               :: coszrs(pcols)
   real(r8), intent(in)               :: solin(pcols)
   real(r8), intent(in)               :: ztodt

   ! SY: for random ensemble averaging
   integer, external :: shuffled_1d 

   ! local variables
   real(r8) :: input(pcols,inputlength)
   real(r8) :: output(pcols,outputlength)
   integer :: i,k,ncol,ixcldice,ixcldliq,ii,kk,idx_trop(pcols),kens
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
      input(:ncol,2*pver+1:3*pver) = state%q(1:ncol,1:pver,ixcldliq) ! state_q0002
      input(:ncol,3*pver+1:4*pver) = state%q(1:ncol,1:pver,ixcldice) ! state_q0003
      input(:ncol,4*pver+1:5*pver) = state%u(1:ncol,1:pver)          ! state_u
      input(:ncol,5*pver+1:6*pver) = state%v(1:ncol,1:pver)          ! state_v
      input(:ncol,6*pver+1       ) = state%ps(1:ncol)                ! state_ps
      input(:ncol,6*pver+2       ) = solin(1:ncol)                   ! pbuf_SOLIN
      input(:ncol,6*pver+3       ) = lhflx(1:ncol)                   ! pbuf_LHFLX
      input(:ncol,6*pver+4       ) = shflx(1:ncol)                   ! pbuf_SHFLX
      input(:ncol,6*pver+5       ) = taux(1:ncol)                    ! pbuf_TAUX
      input(:ncol,6*pver+6       ) = tauy(1:ncol)                    ! pbuf_TAUY
      input(:ncol,6*pver+7       ) = coszrs(1:ncol)                  ! pbuf_COSZRS
      input(:ncol,6*pver+8       ) = cam_in%ALDIF(:ncol)             ! cam_in_ALDIF 
      input(:ncol,6*pver+9       ) = cam_in%ALDIR(:ncol)             ! cam_in_ALDIR 
      input(:ncol,6*pver+10      ) = cam_in%ASDIF(:ncol)             ! cam_in_ASDIF 
      input(:ncol,6*pver+11      ) = cam_in%ASDIR(:ncol)             ! cam_in_ASDIR
      input(:ncol,6*pver+12      ) = cam_in%LWUP(:ncol)              ! cam_in_LWUP
      input(:ncol,6*pver+13      ) = cam_in%ICEFRAC(:ncol)           ! cam_in_ICEFRAC
      input(:ncol,6*pver+14      ) = cam_in%LANDFRAC(:ncol)          ! cam_in_LANDFRAC
      input(:ncol,6*pver+15      ) = cam_in%OCNFRAC(:ncol)           ! cam_in_OCNFRAC
      input(:ncol,6*pver+16      ) = cam_in%SNOWHICE(:ncol)          ! cam_in_SNOWHICE
      input(:ncol,6*pver+17      ) = cam_in%SNOWHLAND(:ncol)         ! cam_in_SNOWHLAND
      input(:ncol,6*pver+18:6*pver+33) = ozone(:ncol,6:21)           ! pbuf_ozone
      input(:ncol,6*pver+34:6*pver+49) = ch4(:ncol,6:21)             ! pbuf_CH4
      input(:ncol,6*pver+50:6*pver+65) = n2o(:ncol,6:21)             ! pbuf_N2O

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

end select

! Tue Jan 24 13:28:43 CST 2023
! Sungduk 
#ifdef CLIMSIMDEBUG
if (masterproc) then ! (logging ncol=1 only)
select case (to_lower(trim(cb_nn_var_combo)))
   case('v2')
      write (iulog,*) 'CLIMSIMDEBUG state_t = ', state%t(1:ncol,1:pver)          ! state_t
      write (iulog,*) 'CLIMSIMDEBUG state_q0001 = ', state%q(1:ncol,1:pver,1)        ! state_q0001
      write (iulog,*) 'CLIMSIMDEBUG state_q0002 = ', state%q(1:ncol,1:pver,ixcldliq) ! state_q0002
      write (iulog,*) 'CLIMSIMDEBUG state_q0003 = ', state%q(1:ncol,1:pver,ixcldice) ! state_q0003
      write (iulog,*) 'CLIMSIMDEBUG state_u = ', state%u(1:ncol,1:pver)          ! state_u
      write (iulog,*) 'CLIMSIMDEBUG state_v = ', state%v(1:ncol,1:pver)          ! state_v
      write (iulog,*) 'CLIMSIMDEBUG state_ps = ', state%ps(1:ncol)                ! state_ps
      write (iulog,*) 'CLIMSIMDEBUG pbuf_SOLIN = ', solin(1:ncol)                   ! pbuf_SOLIN
      write (iulog,*) 'CLIMSIMDEBUG pbuf_LHFLX = ', lhflx(1:ncol)                   ! pbuf_LHFLX
      write (iulog,*) 'CLIMSIMDEBUG pbuf_SHFLX = ', shflx(1:ncol)                   ! pbuf_SHFLX
      write (iulog,*) 'CLIMSIMDEBUG pbuf_TAUX = ', taux(1:ncol)                    ! pbuf_TAUX
      write (iulog,*) 'CLIMSIMDEBUG pbuf_TAUY = ', tauy(1:ncol)                    ! pbuf_TAUY
      write (iulog,*) 'CLIMSIMDEBUG pbuf_COSZRS = ', coszrs(1:ncol)                  ! pbuf_COSZRS
      write (iulog,*) 'CLIMSIMDEBUG cam_in_ALDIF = ', cam_in%ALDIF(:ncol)             ! cam_in_ALDIF
      write (iulog,*) 'CLIMSIMDEBUG cam_in_ALDIR = ', cam_in%ALDIR(:ncol)             ! cam_in_ALDIR
      write (iulog,*) 'CLIMSIMDEBUG cam_in_ASDIF = ', cam_in%ASDIF(:ncol)             ! cam_in_ASDIF
      write (iulog,*) 'CLIMSIMDEBUG cam_in_ASDIR = ', cam_in%ASDIR(:ncol)             ! cam_in_ASDIR
      write (iulog,*) 'CLIMSIMDEBUG cam_in_LWUP = ', cam_in%LWUP(:ncol)              ! cam_in_LWUP
      write (iulog,*) 'CLIMSIMDEBUG cam_in_ICEFRAC = ', cam_in%ICEFRAC(:ncol)           ! cam_in_ICEFRAC
      write (iulog,*) 'CLIMSIMDEBUG cam_in_LANDFRAC = ', cam_in%LANDFRAC(:ncol)          ! cam_in_LANDFRAC
      write (iulog,*) 'CLIMSIMDEBUG cam_in_OCNFRAC = ', cam_in%OCNFRAC(:ncol)           ! cam_in_OCNFRAC
      write (iulog,*) 'CLIMSIMDEBUG cam_in_SNOWHICE = ', cam_in%SNOWHICE(:ncol)          ! cam_in_SNOWHICE
      write (iulog,*) 'CLIMSIMDEBUG cam_in_SNOWHLAND = ', cam_in%SNOWHLAND(:ncol)         ! cam_in_SNOWHLAND
      write (iulog,*) 'CLIMSIMDEBUG pbuf_ozone = ', ozone(:ncol,6:21)           ! pbuf_ozone
      write (iulog,*) 'CLIMSIMDEBUG pbuf_CH4 = ', ch4(:ncol,6:21)             ! pbuf_CH4
      write (iulog,*) 'CLIMSIMDEBUG pbuf_N2O = ', n2o(:ncol,6:21)             ! pbuf_N2O
      if (input_rh) then ! relative humidity conversion for input
         write (iulog,*) 'CLIMSIMDEBUG RH = ', input(:ncol,1*pver+1:2*pver) ! relhum 
      end if
end select
end if
#endif 

#ifdef CLIMSIMDEBUG
      if (masterproc) then
        write (iulog,*) 'CLIMSIMDEBUG input pre norm=',input(1,:)
      endif
#endif

    ! 2. Normalize input
    do k=1,inputlength
      input(:ncol,k) = (input(:ncol,k) - inp_sub(k))/inp_div(k)
    end do
#ifdef CLIMSIMDEBUG
      if (masterproc) then
        write (iulog,*) 'CLIMSIMDEBUG input post norm=',input(1,:)
      endif
#endif

    do i=1,ncol
      if (cb_do_ensemble) then
        output(i,:) = 0.
        !! Random ensemble averaging
        if (cb_do_random_ensemble) then
          ens_ind_shuffled = shuffle_1d(ens_ind_shuffled) ! randomly shuffle ens indices
          do kens = 1,cb_random_ens_size
            output(i,:) = output(i,:) +  (1._r8/cb_random_ens_size) * climsim_net(ens_ind_shuffled(kens)) % output(input(i,:))
          enddo
#ifdef CLIMSIMDEBUG
          if (masterproc .and. i.eq.1) then
            write (iulog,*) 'CLIMSIMDEBUG  random ensemble model IDs = ',ens_ind_shuffled(1:cb_random_ens_size)
          endif
#endif
        !! All ensemble averaging
        else
          do kens = 1,cb_ens_size
            output(i,:) = output(i,:) +  (1._r8/cb_ens_size) * climsim_net(kens) % output(input(i,:))
          enddo
        endif
      !! Using a single model
      else ! cb_do_ensemble
        output(i,:) = climsim_net(1) % output(input(i,:))
      endif
    end do
#ifdef CLIMSIMDEBUG
      if (masterproc) then
        write (iulog,*) 'CLIMSIMDEBUG output = ',output(1,:)
      endif
#endif

   ! Manually applying ReLU activation for positive-definite variables
   ! [TODO] for ensemble, ReLU should be moved before ens-averaging
   do i=1,ncol
     ! ReLU for the last 8 variables (true for 'v1' and 'v2')
     do k=outputlength-7,outputlength
       output(i,k) = max(output(i,k), 0.)
     end do
     ! tiny flwds
     k=4*pver+3
     output(i,k) = max(output(i,k), tiny(output(i,k))) ! flwds
                                                       ! preventing flwds==0 error
     ! zero out surface solar fluxes when local time is at night
     if (coszrs(i) .le. 0.) then
       output(i,6*pver+1) = 0. ! netsw
       output(i,6*pver+5) = 0. ! sols
       output(i,6*pver+6) = 0. ! soll
       output(i,6*pver+7) = 0. ! solsd
       output(i,6*pver+8) = 0. ! solld
     endif
   end do
#ifdef CLIMSIMDEBUG
      if (masterproc) then
        write (iulog,*) 'CLIMSIMDEBUG output after ReLU = ',output(1,:)
      endif
#endif

   ! output normalization (un-weighting, really).
   do i=1,ncol
     do k=1,outputlength
      output(i,k) = output(i,k) / out_scale(k)
     end do
   end do
#ifdef CLIMSIMDEBUG
      if (masterproc) then
        write (iulog,*) 'CLIMSIMDEBUG output post scale = ',output(1,:)
      endif
#endif

! ---------- 1. NN output to atmosphere forcing --------
! ['TBCTEND', 'QBCTEND','CLDLIQBCTEND','CLDICEBCTEND']
   s_bctend (:ncol,1:pver) = output(1:ncol,0*pver+1:1*pver)*cpair ! K/s --> J/kg/s (ptend expects that)
   q_bctend (:ncol,1:pver) = output(1:ncol,1*pver+1:2*pver)       ! kg/kg/s
   qc_bctend(:ncol,1:pver) = output(1:ncol,2*pver+1:3*pver)       ! kg/kg/s
   qi_bctend(:ncol,1:pver) = output(1:ncol,3*pver+1:4*pver)       ! kg/kg/s
   u_bctend (:ncol,1:pver) = output(1:ncol,4*pver+1:5*pver)       ! m/s/s
   v_bctend (:ncol,1:pver) = output(1:ncol,5*pver+1:6*pver)       ! m/s/s

! deny any moisture activity in the stratosphere:
   do i=1,ncol
     call detect_tropopause(state%t(i,:),state%exner(i,:),state%zm(i,:),state%pmid(i,:),idx_trop(i))
     q_bctend (i,1:idx_trop(i)) = 0.
     qc_bctend(i,1:idx_trop(i)) = 0.
     qi_bctend(i,1:idx_trop(i)) = 0.
   end do
   call outfld('TROP_IND', idx_trop(:ncol)*1._r8, ncol, state%lchnk)

! -- atmos positivity constraints ---- 
   if (do_constraints) then
   do i=1,ncol
     do k=1,pver
! deny activity in the ice phase where it is above freezing.
       if (state%t(i,k) .gt. 273.16) then
          qi_bctend(i,k) = 0.
! deny activitiy in the water phase where it is below freezing.
! (253.16K: the lowest threshold temperature for supercooled cloud water to form)
       elseif (state%t(i,k) .lt. 253.16) then
          qc_bctend(i,k) = 0.
       end if
!eliminate all activity in the water phase on top 10 levels:

! energy positivity:
       safter = state%s(i,k) + s_bctend(i,k)*ztodt ! predicted DSE after NN tendency
       if (safter .lt. 0.) then ! can only happen when bctend < 0...
         s_bctend(i,k) = s_bctend(i,k) + abs(safter)/ztodt ! in which case reduce cooling rate
         write (iulog,*) 'HEY CLIMSIM made a negative absolute temperature, corrected but BEWARE!!!'
         write (iulog,*) '' ! [TODO] printout lat/lon and error magnitude
       endif

 ! vapor positivity:
       qafter = state%q(i,k,1) + q_bctend(i,k)*ztodt ! predicted vapor after NN tendency
       if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
         q_bctend(i,k) = q_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
         write (iulog,*) 'HEY CLIMSIM made a negative absolute q, corrected but BEWARE!!!'
       endif

 ! liquid positivity:
       qafter = state%q(i,k,ixcldliq) + qc_bctend(i,k)*ztodt ! predicted liquid after NN tendency
       if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
         qc_bctend(i,k) = qc_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
         write (iulog,*) 'HEY CLIMSIM made a negative absolute qc, corrected but BEWARE!!!'
       endif
! ice positivity:
       qafter = state%q(i,k,ixcldice) + qi_bctend(i,k)*ztodt ! predicted ice after NN tendency
       if (qafter .lt. 0.) then ! can only happen when qbctend < 0...
         qi_bctend(i,k) = qi_bctend(i,k) + abs(qafter)/ztodt ! in which case reduce drying rate
         write (iulog,*) 'HEY CLIMSIM made a negative absolute qi, corrected but BEWARE!!!'
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
!!!          (#CLIMSIM_OCN_ONLY)
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
#ifdef CLIMSIM_OCN_ONLY 
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
#ifdef CLIMSIM_OCN_ONLY
     end if
#endif
   end do 

end subroutine neural_net

  subroutine init_neural_net()

    implicit none

    integer :: i, k

    allocate(inp_sub (inputlength))
    allocate(inp_div (inputlength))
    allocate(out_scale (outputlength))
    
    ! ens-mean inference
    if (cb_do_ensemble) then
       write (iulog,*) 'CLIMSIM: Ensemble is turned on with Ensemble size  ', cb_ens_size
       allocate(climsim_net (cb_ens_size))
       do i = 1,cb_ens_size
          call climsim_net(i) %load(cb_ens_fkb_model_list(i))
          write (iulog,*) 'CLIMSIM: Ensemble fkb model (', i, ') : ', trim(cb_ens_fkb_model_list(i))
       enddo

       ! random ensemble
       if (cb_random_ens_size .ge. 1) then
          write (iulog,*) 'CLIMSIM: Random ensemble averaging with N = ', cb_random_ens_size
          if (cb_random_ens_size .le. cb_ens_size) then
             allocate(ens_ind_shuffled(cb_ens_size))
             ens_ind_shuffled = (/ (k, k=1, cb_ens_size) /)
             cb_do_random_ensemble = .true.
          else
             call endrun("init_neural_net error: cb_random_ens_size should be less than or equal to cb_ens_size")
          endif

       endif

    ! single model inference
    else
       allocate(climsim_net (1))
       call climsim_net(1) %load(cb_fkb_model)
       if (masterproc) then
          write (iulog,*) 'CLIMSIM: loaded network from txt file, ', trim(cb_fkb_model)
       endif
    endif

    open (unit=555,file=cb_inp_sub,status='old',action='read')
    read(555,*) inp_sub(:)
    close (555)
    if (masterproc) then
       write (iulog,*) 'CLIMSIM: loaded input subtraction factors from: ', trim(cb_inp_sub)
    endif

    open (unit=555,file=cb_inp_div,status='old',action='read')
    read(555,*) inp_div(:)
    close (555)
    if (masterproc) then
       write (iulog,*) 'CLIMSIM: loaded input division factors from: ', trim(cb_inp_div)
    endif

    open (unit=555,file=cb_out_scale,status='old',action='read')
    read(555,*) out_scale(:)
    close (555)
    if (masterproc) then
       write (iulog,*) 'CLIMSIM: loaded output scale factors from: ', trim(cb_out_scale)
    endif

#ifdef CLIMSIMDEBUG
    if (masterproc) then
       write (iulog,*) 'CLIMSIMDEBUG read input norm inp_sub=', inp_sub(:)
       write (iulog,*) 'CLIMSIMDEBUG read input norm inp_div=', inp_div(:)       
       write (iulog,*) 'CLIMSIMDEBUG read output norm out_scale=', out_scale(:)       
    endif
#endif

  ! add diagnostic output fileds
  call addfld ('TROP_IND',horiz_only,   'A', '1', 'lev index for tropopause')

  end subroutine init_neural_net

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
      tom_eice = 100.*(cice(4) + max(cice(2),dt)*(cice(5)+max(cice(3),dt)*cice(6))) 
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
     if (pmid(k) .le. 40000. .and. dthetadz > 10.) then
       klev_crmtop = k
     endif
   end do
  end subroutine detect_tropopause

  ! Read namelist variables.
  subroutine climsim_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr, f
      character(len=*), parameter :: subname = 'climsim_readnl'
      
      namelist /climsim_nl/ inputlength, outputlength, input_rh, &
                           cb_fkb_model, &
                           cb_inp_sub, cb_inp_div, cb_out_scale, &
                           cb_partial_coupling, cb_partial_coupling_vars,&
                           cb_use_input_prectm1, &
                           cb_do_ensemble, cb_ens_size, cb_ens_fkb_model_list, &
                           cb_random_ens_size, &
                           cb_nn_var_combo

      ! Initialize 'cb_partial_coupling_vars'
      do f = 1, pflds
        cb_partial_coupling_vars(f) = ' '
      end do

      ! Initialize 'cb_ens_fkb_model_list'
      allocate(cb_ens_fkb_model_list(max_nn_ens))
      do f = 1, max_nn_ens
        cb_ens_fkb_model_list(f) = ' '
      end do

      ierr = 0
      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'climsim_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, climsim_nl, iostat=ierr)
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
      call mpibcast(cb_fkb_model, len(cb_fkb_model), mpichar, 0, mpicom)
      call mpibcast(cb_inp_sub,   len(cb_inp_sub),   mpichar, 0, mpicom)
      call mpibcast(cb_inp_div,   len(cb_inp_div),   mpichar, 0, mpicom)
      call mpibcast(cb_out_scale, len(cb_out_scale), mpichar, 0, mpicom)
      call mpibcast(cb_partial_coupling, 1,          mpilog,  0, mpicom)
      call mpibcast(cb_partial_coupling_vars, len(cb_partial_coupling_vars(1))*pflds, mpichar, 0, mpicom)
      call mpibcast(cb_do_ensemble, 1,               mpilog,  0, mpicom)
      call mpibcast(cb_ens_size,    1,               mpiint,  0, mpicom)
      call mpibcast(cb_ens_fkb_model_list,    len(cb_ens_fkb_model_list(1))*max_nn_ens, mpichar, 0, mpicom)
      call mpibcast(cb_random_ens_size,    1,        mpiint,  0, mpicom)
      call mpibcast(cb_nn_var_combo, len(cb_nn_var_combo), mpichar,  0, mpicom)
      ! [TODO] check ierr for each mpibcast call
      ! if (ierr /= 0) then
      !    call endrun(subname // ':: ERROR broadcasting namelist variable cb_partial_coupling_vars')
      ! end if
#endif

   end subroutine climsim_readnl

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

end module climsim
