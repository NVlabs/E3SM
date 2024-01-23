#ifdef MMF_ML_TRAINING
module ml_training

   ! TODO:
   ! - figure out how to handle prescribed aerosol, ozone, etc. (see phys_timestep_init)
   ! - allow data to be written at specific interval (i.e. not every time step)

   use shr_kind_mod,       only: r8 => shr_kind_r8
   use spmd_utils,         only: masterproc
   use constituents,       only: pcnst ! for q and cflx
   use seq_drydep_mod,     only: n_drydep ! for depvel
   use shr_megan_mod,      only: shr_megan_mechcomps_n ! for meganflx
   use ppgrid,             only: pver, pverp, pcols, begchunk, endchunk
   use cam_abortutils,     only: endrun
   use cam_history_support,only: fillvalue
   use camsrfexch,         only: cam_in_t, cam_out_t
   use cam_logfile,        only: iulog
   use shr_sys_mod,        only: shr_sys_flush
   use physics_types,      only: physics_state, physics_tend
   use physics_buffer,     only: physics_buffer_desc
   use pio,                only: file_desc_t, io_desc_t, var_desc_t, &
                                 pio_double, pio_int, pio_noerr, &
                                 pio_seterrorhandling, pio_bcast_error, &
                                 pio_def_var, pio_def_dim, pio_enddef, &
                                 pio_put_var, pio_write_darray, pio_closefile
   use perf_mod,           only: t_startf, t_stopf

   implicit none
   private
   save
   
   ! Public interfaces
   public :: write_ml_training

   ! internal private routines
   private :: get_ml_filename

   integer, parameter :: filename_len = 256

CONTAINS
   !------------------------------------------------------------------------------------------------
   function get_ml_filename(fspec_in,yr,mn,dy,sec) result(fname)
      ! return filename based on date/time info and fspec_in specifiers
      use seq_timemgr_mod, only: seq_timemgr_EClockGetData
      use filenames,       only: interpret_filename_spec
      character(len=4), intent(in) :: fspec_in
      integer,          intent(in) :: yr,mn,dy,sec   ! current year, month, day, and time of day (sec)
      character(len=filename_len) :: fspec_loc
      character(len=filename_len) :: fname  ! full file name to return
      ! (%c=caseid, $y=year, $m=month, $d=day, $s=sec in day, %t=number)
      fspec_loc = '%c.eam.'//trim(fspec_in)//'.%y-%m-%d-%s.nc'
      fname = interpret_filename_spec( fspec_loc, yr_spec=yr, mon_spec=mn, day_spec=dy, sec_spec=sec )
   end function get_ml_filename
   !------------------------------------------------------------------------------------------------
   subroutine write_ml_training( pbuf2d, phys_state, phys_tend, cam_in, cam_out, yr, mn, dy, sec, mode )
      use phys_grid,           only: phys_decomp
      use physics_buffer,      only: pbuf_init_restart_alt, pbuf_write_restart_alt
      use time_manager,        only: timemgr_init_restart, timemgr_write_restart
      use chemistry,           only: chem_init_restart, chem_write_restart
      use cam_grid_support,    only: cam_grid_id, cam_grid_header_info_t
      use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
      use cam_grid_support,    only: cam_grid_write_attr, cam_grid_write_var
      use cam_pio_utils,       only: cam_pio_createfile, cam_pio_def_dim, cam_pio_closefile
      use phys_control,        only: phys_getopts
      use hycoef,              only: init_restart_hycoef
      use dyn_grid,            only: get_horiz_grid_d
      !-------------------------------------------------------------------------
      ! Input arguments
      type(physics_buffer_desc), pointer    :: pbuf2d(:,:)
      type(physics_state),       intent(in) :: phys_state(begchunk:endchunk)
      type(physics_tend ),       intent(in) :: phys_tend (begchunk:endchunk)
      type(cam_in_t),            intent(in) :: cam_in    (begchunk:endchunk)
      type(cam_out_t),           intent(in) :: cam_out   (begchunk:endchunk)
      integer,                   intent(in) :: yr   ! Simulation year
      integer,                   intent(in) :: mn   ! Simulation month
      integer,                   intent(in) :: dy   ! Simulation day
      integer,                   intent(in) :: sec  ! Seconds into current simulation day
      integer,                   intent(in) :: mode ! used to select input / output modes
      !-------------------------------------------------------------------------
      ! Local workspace
      type(file_desc_t)            :: file
      integer                      :: pver_id, pverp_id, pcnst_id
      integer                      :: ncol_dimid
      integer                      :: grid_id
      type(cam_grid_header_info_t) :: header_info ! A structure to hold the horz dims and coord info

      integer :: ncol(begchunk:endchunk)                ! ncol value per chunk
      type(io_desc_t), pointer :: iodesc2d              ! pcols * (endchunk-begchunk+1)
      type(io_desc_t), pointer :: iodesc3d              ! pcols * pver  * (endchunk-begchunk+1) 
      type(io_desc_t), pointer :: iodesc3dp             ! pcols * pverp * (endchunk-begchunk+1)
      real(r8):: tmp2D(pcols, begchunk:endchunk)          ! temp variable for derived type data
      real(r8):: tmp3D(pcols, pver, begchunk:endchunk)    ! temp variable for derived type data
      real(r8):: tmp3Dp(pcols, pverp, begchunk:endchunk)  ! temp variable for derived type data

      integer :: ierr, i, m
      integer :: physgrid
      integer :: gdims(3)
      integer :: nhdims

      integer, parameter, dimension(1) :: dimids_hrz = (/1/)     ! horz dims only
      integer, parameter, dimension(2) :: dimids_3D1 = (/1,2/)   ! horz + pver
      integer, parameter, dimension(2) :: dimids_3D2 = (/1,3/)   ! horz + pverp

      character(len=8)     :: num      ! used for writing numeric charaters (i.e. constituent index)
      character(len=4)     :: fspec    ! string used after ".eam." in file name 
      logical              :: add_pbuf
      logical              :: add_phys_state
      logical              :: add_phys_tend
      logical              :: add_cam_in
      logical              :: add_cam_out

      ! file variable descriptions

      ! cam_out
      ! [0] tbot(:)     ! bot level temperature
      ! [0] zbot(:)     ! bot level height above surface
      ! [0] ubot(:)     ! bot level u wind
      ! [0] vbot(:)     ! bot level v wind
      ! [0] qbot(:,:)   ! bot level specific humidity
      ! [0] pbot(:)     ! bot level pressure
      ! [0] rho(:)      ! bot level density
      ! [1] netsw(:)    !
      ! [1] flwds(:)    !
      ! [1] precsc(:)   !
      ! [0] precsl(:)   !
      ! [1] precc(:)    !
      ! [0] precl(:)    !
      ! [1] soll(:)     !
      ! [1] sols(:)     !
      ! [1] solld(:)    !
      ! [1] solsd(:)    !
      ! [0] thbot(:)    !
      ! [0] co2prog(:)  ! prognostic co2
      ! [0] co2diag(:)  ! diagnostic co2
      ! [0] psl(:)      ! sea level pressure
      ! [0] bcphiwet(:) ! wet deposition of hydrophilic black carbon
      ! [0] bcphidry(:) ! dry deposition of hydrophilic black carbon
      ! [0] bcphodry(:) ! dry deposition of hydrophobic black carbon
      ! [0] ocphiwet(:) ! wet deposition of hydrophilic organic carbon
      ! [0] ocphidry(:) ! dry deposition of hydrophilic organic carbon
      ! [0] ocphodry(:) ! dry deposition of hydrophobic organic carbon
      ! [0] dstwet1(:)  ! wet deposition of dust (bin1)
      ! [0] dstdry1(:)  ! dry deposition of dust (bin1)
      ! [0] dstwet2(:)  ! wet deposition of dust (bin2)
      ! [0] dstdry2(:)  ! dry deposition of dust (bin2)
      ! [0] dstwet3(:)  ! wet deposition of dust (bin3)
      ! [0] dstdry3(:)  ! dry deposition of dust (bin3)
      ! [0] dstwet4(:)  ! wet deposition of dust (bin4)
      ! [0] dstdry4(:)  ! dry deposition of dust (bin4)
      ! [0] wsresp(:)   ! first-order response of low-level wind to surface fluxes
      ! [0] tau_est(:)  ! stress estimated to be in equilibrium with ubot/vbot
      ! [0] ugust(:)    ! gustiness value
      ! [0] uovern(:)       ! ratio of wind speed/brunt vaisalla frequency
      ! type(var_desc_t)     :: desc_trefmxav
      ! type(var_desc_t)     :: desc_trefmnav
      ! type(var_desc_t)     :: desc_tbot
      ! type(var_desc_t)     :: desc_zbot
      ! type(var_desc_t)     :: desc_ubot
      ! type(var_desc_t)     :: desc_vbot
      ! type(var_desc_t)     :: desc_qbot(pcnst)
      ! type(var_desc_t)     :: desc_pbot
      ! type(var_desc_t)     :: desc_rho
      type(var_desc_t)     :: desc_netsw
      type(var_desc_t)     :: desc_flwds
      type(var_desc_t)     :: desc_precsc
      ! type(var_desc_t)     :: desc_precsl
      type(var_desc_t)     :: desc_precc
      ! type(var_desc_t)     :: desc_precl
      type(var_desc_t)     :: desc_solld
      type(var_desc_t)     :: desc_sols
      type(var_desc_t)     :: desc_soll
      type(var_desc_t)     :: desc_solsd
      ! type(var_desc_t)     :: desc_thbot
      ! type(var_desc_t)     :: desc_psl
      ! type(var_desc_t)     :: desc_wsresp
      ! type(var_desc_t)     :: desc_tau_est
      ! type(var_desc_t)     :: desc_ugust
      ! type(var_desc_t)     :: desc_uovern
      ! type(var_desc_t)     :: desc_co2prog
      ! type(var_desc_t)     :: desc_co2diag
      ! type(var_desc_t)     :: desc_bcphidry ! surface deposition flux are not used when atm_dep_flux=.false. in atm_in
      ! type(var_desc_t)     :: desc_bcphodry
      ! type(var_desc_t)     :: desc_ocphidry
      ! type(var_desc_t)     :: desc_ocphodry
      ! type(var_desc_t)     :: desc_dstdry1
      ! type(var_desc_t)     :: desc_dstdry2
      ! type(var_desc_t)     :: desc_dstdry3
      ! type(var_desc_t)     :: desc_dstdry4
      ! type(var_desc_t)     :: desc_bcphiwet
      ! type(var_desc_t)     :: desc_ocphiwet
      ! type(var_desc_t)     :: desc_dstwet1
      ! type(var_desc_t)     :: desc_dstwet2
      ! type(var_desc_t)     :: desc_dstwet3
      ! type(var_desc_t)     :: desc_dstwet4

      ! physics_state
      ! [1] ps(:)        ! surface pressure
      ! [0] psdry(:)     ! dry surface pressure
      ! [0] phis(:)      ! surface geopotential
      ! [0] ulat(:)      ! unique latitudes  (radians)
      ! [0] ulon(:)      ! unique longitudes (radians)
      ! [1] t(:,:)       ! temperature (K)
      ! [1] u(:,:)       ! zonal wind (m/s)
      ! [1] v(:,:)       ! meridional wind (m/s)
      ! [0] s(:,:)       ! dry static energy
      ! [0] omega(:,:)   ! vertical pressure velocity (Pa/s)
      ! [1] pmid(:,:)    ! midpoint pressure (Pa)
      ! [0] pmiddry(:,:) ! midpoint pressure dry (Pa)
      ! [0] pdel(:,:)    ! layer thickness (Pa)
      ! [0] pdeldry(:,:) ! layer thickness dry (Pa)
      ! [0] rpdel(:,:)   ! reciprocal of layer thickness (Pa)
      ! [0] rpdeldry(:,:)! recipricol layer thickness dry (Pa)
      ! [0] lnpmid(:,:)  ! ln(pmid)
      ! [0] lnpmiddry(:,:)! log midpoint pressure dry (Pa)
      ! [0] exner(:,:)   ! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
      ! [0] zm(:,:)      ! geopotential height above surface at midpoints (m)
      ! [1] q (:,:,:)    ! constituent mixing ratio (kg/kg moist or dry air depending on type)
      !                  ! Only print out q000{1,2,3} that accounts for Q, CLDLIQ, and CLDICE.
      !
      ! (Diagnostic, grid, carbon-flux variables are omitted:
      !  pint, pintdry, lnpint, lnpintdry, zi,
      !  te_ini, te_cur, tw_ini, tw_cur, tc_curr, tc_init, tc_mnst, tc_prev,
      !  c_flux_sfc, c_mflx_air, c_mflx_sff, c_mflx_lnd, c_mflx_ocn,
      !  c_iflx_sfc, c_iflx_air, c_iflx_sff, c_iflx_lnd, c_iflx_ocn
      ! )

      type(var_desc_t)     :: state_desc_ps
      ! type(var_desc_t)     :: state_desc_psdry
      ! type(var_desc_t)     :: state_desc_phis
      type(var_desc_t)     :: state_desc_t
      type(var_desc_t)     :: state_desc_u
      type(var_desc_t)     :: state_desc_v
      type(var_desc_t)     :: state_desc_s
      ! type(var_desc_t)     :: state_desc_omega
      type(var_desc_t)     :: state_desc_pmid
      ! type(var_desc_t)     :: state_desc_pmiddry
      ! type(var_desc_t)     :: state_desc_pdel
      ! type(var_desc_t)     :: state_desc_pdeldry
      ! type(var_desc_t)     :: state_desc_rpdel
      ! type(var_desc_t)     :: state_desc_rpdeldry
      ! type(var_desc_t)     :: state_desc_lnpmid
      ! type(var_desc_t)     :: state_desc_lnpmiddry
      ! type(var_desc_t)     :: state_desc_exner
      ! type(var_desc_t)     :: state_desc_zm
      type(var_desc_t)     :: state_desc_q(3) !state_desc_q(pcnst)
      ! type(var_desc_t)     :: state_desc_pint
      ! type(var_desc_t)     :: state_desc_pintdry
      ! type(var_desc_t)     :: state_desc_lnpint
      ! type(var_desc_t)     :: state_desc_lnpintdry
      ! type(var_desc_t)     :: state_desc_zi

      ! physics_tend
      ! [0] dtdt(:,:)
      ! [0] dudt(:,:)
      ! [0] dvdt(:,:)
      ! [0] flx_net(:)
      ! [0] te_tnd(:) ! cumulative boundary flux of total energy
      ! [0] tw_tnd(:) ! cumulative boundary flux of total water
      ! type(var_desc_t)     :: tend_desc_dtdt
      ! type(var_desc_t)     :: tend_desc_dudt
      ! type(var_desc_t)     :: tend_desc_dvdt
      ! type(var_desc_t)     :: tend_desc_flx_net
      ! type(var_desc_t)     :: tend_desc_te_tnd
      ! type(var_desc_t)     :: tend_desc_tw_tnd
      
      ! cam_in
      ! [1] asdir(:)      ! albedo: shortwave, direct
      ! [1] asdif(:)      ! albedo: shortwave, diffuse
      ! [1] aldir(:)      ! albedo: longwave, direct
      ! [1] aldif(:)      ! albedo: longwave, diffuse
      ! [1] lwup(:)       ! longwave up radiative flux
      ! [0] lhf(:)        ! latent heat flux
      ! [0] shf(:)        ! sensible heat flux
      ! [0] h2otemp(:)    ! water temperature heat flux from ocean
      ! [0] wsx(:)        ! surface u-stress (N)
      ! [0] wsy(:)        ! surface v-stress (N)
      ! [0] tref(:)       ! ref height surface air temp
      ! [0] qref(:)       ! ref height specific humidity
      ! [0] u10(:)        ! 10m wind speed
      ! [0] ts(:)         ! merged surface temp
      ! [0] sst(:)        ! sea surface temp
      ! [1] snowhland(:)  ! snow depth (liquid water equivalent) over land
      ! [1] snowhice(:)   ! snow depth over ice
      ! [0] fco2_lnd(:)   ! co2 flux from lnd
      ! [0] fco2_ocn(:)   ! co2 flux from ocn
      ! [0] fdms(:)       ! dms flux
      ! [1] landfrac(:)   ! land area fraction
      ! [1] icefrac(:)    ! sea-ice areal fraction
      ! [1] ocnfrac(:)    ! ocean areal fraction
      ! [0] ram1       !aerodynamical resistance (s/m) (pcols)
      ! [0] fv         !friction velocity (m/s) (pcols)
      ! [0] soilw      !volumetric soil water (m3/m3)
      ! [0] cflx(:,:)     ! constituent flux (emissions), dims: [:,pcnst]
      ! [0] ustar(:)      ! atm/ocn saved version of ustar
      ! [0] re(:)         ! atm/ocn saved version of re
      ! [0] ssq(:)        ! atm/ocn saved version of ssq
      ! [0] depvel(:,:)   ! deposition velocities, dims: [:,n_drydep]
      ! [0] dstflx(:,:)   ! dust fluxes, dims: [:,4] ! 4 bins from surface model
      ! [0] meganflx(:,:) ! MEGAN fluxes, dims: [:,shr_megan_mechcomps_n]
      type(var_desc_t)     :: desc_asdir
      type(var_desc_t)     :: desc_asdif
      type(var_desc_t)     :: desc_aldir
      type(var_desc_t)     :: desc_aldif
      type(var_desc_t)     :: desc_lwup
      ! type(var_desc_t)     :: desc_lhf
      ! type(var_desc_t)     :: desc_shf
      ! type(var_desc_t)     :: desc_h2otemp
      ! type(var_desc_t)     :: desc_wsx
      ! type(var_desc_t)     :: desc_wsy
      ! type(var_desc_t)     :: desc_tref
      ! type(var_desc_t)     :: desc_qref
      ! type(var_desc_t)     :: desc_u10
      ! type(var_desc_t)     :: desc_ts
      ! type(var_desc_t)     :: desc_sst
      type(var_desc_t)     :: desc_snowhland
      type(var_desc_t)     :: desc_snowhice
      ! type(var_desc_t)     :: desc_fco2_lnd
      ! type(var_desc_t)     :: desc_fco2_ocn
      ! type(var_desc_t)     :: desc_fdms
      type(var_desc_t)     :: desc_landfrac
      type(var_desc_t)     :: desc_icefrac
      type(var_desc_t)     :: desc_ocnfrac
      ! type(var_desc_t)     :: desc_ram1
      ! type(var_desc_t)     :: desc_fv
      ! type(var_desc_t)     :: desc_soilw ! MOZART
      ! type(var_desc_t)     :: desc_cflx(pcnst)
      ! type(var_desc_t)     :: desc_ustar
      ! type(var_desc_t)     :: desc_re
      ! type(var_desc_t)     :: desc_ssq
      ! type(var_desc_t)     :: desc_depvel(n_drydep)
      ! type(var_desc_t)     :: desc_dstflx(4)
      ! type(var_desc_t)     :: desc_meganflx(shr_megan_mechcomps_n)
      !-------------------------------------------------------------------------
      ! Initialize stuff
      !-------------------------------------------------------------------------
      add_pbuf       = .false.
      add_phys_state = .false.
      add_phys_tend  = .false.
      add_cam_in     = .false.
      add_cam_out    = .false.
      if (mode==1) then
         fspec = 'mli'
         add_pbuf        = .true.
         add_phys_state  = .true.
         add_cam_in      = .true.
      end if
      if (mode==2) then
         fspec = 'mlo'
         add_phys_state = .true.
         add_cam_out    = .true.
      end if

      do i=begchunk,endchunk
         ncol(i) = phys_state(i)%ncol
      end do

      grid_id = cam_grid_id('physgrid')
      call cam_grid_dimensions(grid_id, gdims(1:2), nhdims)
      if (nhdims==1) then
         call cam_grid_get_decomp(grid_id, (/pcols,endchunk-begchunk+1/), (/gdims(1)/),      pio_double, iodesc2d)
         call cam_grid_get_decomp(grid_id, (/pcols,pver,endchunk-begchunk+1/), (/gdims(1),pver/), pio_double, iodesc3d)
         call cam_grid_get_decomp(grid_id, (/pcols,pverp,endchunk-begchunk+1/), (/gdims(1),pverp/), pio_double, iodesc3dp)
      end if
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Initialize file and define variables
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      call cam_pio_createfile(file, trim(get_ml_filename(fspec,yr,mn,dy,sec)))      
      call cam_grid_write_attr(file, grid_id, header_info)

      call cam_pio_def_dim(file, 'lev',   pver,  pver_id,  existOK=.true.)
      call cam_pio_def_dim(file, 'ilev',  pverp, pverp_id, existOK=.true.)
      !call cam_pio_def_dim(file, 'pcnst', pcnst, pcnst_id, existOK=.true.)
      call cam_pio_def_dim(file, 'pcnst', 3, pcnst_id, existOK=.true.)

      call timemgr_init_restart(File)

      if (add_pbuf) call pbuf_init_restart_alt(file, pbuf2d)

      ! data for prognostic chemistry (does nothing for "chem none" option)
      call chem_init_restart(file)
      
      !-------------------------------------------------------------------------
      ! define physics state variables
      if (add_phys_state) then
         if (mode==1) then
           ierr = pio_def_var(file, 'state_ps',        pio_double, dimids_hrz, state_desc_ps)
         end if
         ! ierr = pio_def_var(file, 'state_psdry',     pio_double, dimids_hrz, state_desc_psdry)
         ! ierr = pio_def_var(file, 'state_phis',        pio_double, dimids_hrz, state_desc_phis)
         do m=1,3 !m=1,pcnst
            write(num,'(i4.4)') m
            ierr = pio_def_var(file, 'state_q'//num, pio_double, dimids_3D1, state_desc_q(m))
         end do
         ierr = pio_def_var(file, 'state_t',         pio_double, dimids_3D1, state_desc_t)
         !! SY: future ref
         !    call cam_pio_def_var(File, trim(hist_coords(mdimind)%name), dtype,    &
         !                         (/dimid/), vardesc, existOK=.false.)
         !   ! long_name
         !   ierr=pio_put_att(File, vardesc, 'long_name', trim(hist_coords(mdimind)%long_name))
         !   call cam_pio_handle_error(ierr, 'Error writing "long_name" attr in write_hist_coord_attr')
         !! SY
         ierr = pio_def_var(file, 'state_u',         pio_double, dimids_3D1, state_desc_u)
         ierr = pio_def_var(file, 'state_v',         pio_double, dimids_3D1, state_desc_v)
         ! ierr = pio_def_var(file, 'state_s',         pio_double, dimids_3D1, state_desc_s)
         ! ierr = pio_def_var(file, 'state_omega',     pio_double, dimids_3D1, state_desc_omega)
         if (mode==1) then
           ierr = pio_def_var(file, 'state_pmid',      pio_double, dimids_3D1, state_desc_pmid)
         end if
         ! ierr = pio_def_var(file, 'state_pmiddry',   pio_double, dimids_3D1, state_desc_pmiddry)
         ! ierr = pio_def_var(file, 'state_pdel',      pio_double, dimids_3D1, state_desc_pdel)
         ! ierr = pio_def_var(file, 'state_pdeldry',   pio_double, dimids_3D1, state_desc_pdeldry)
         ! ierr = pio_def_var(file, 'state_rpdel',     pio_double, dimids_3D1, state_desc_rpdel)
         ! ierr = pio_def_var(file, 'state_rpdeldry',  pio_double, dimids_3D1, state_desc_rpdeldry)
         ! ierr = pio_def_var(file, 'state_lnpmid',    pio_double, dimids_3D1, state_desc_lnpmid)
         ! ierr = pio_def_var(file, 'state_lnpmiddry', pio_double, dimids_3D1, state_desc_lnpmiddry)
         ! ierr = pio_def_var(file, 'state_exner',     pio_double, dimids_3D1, state_desc_exner)
         ! ierr = pio_def_var(file, 'state_zm',        pio_double, dimids_3D1, state_desc_zm)
         ! ierr = pio_def_var(file, 'state_pint',      pio_double, dimids_3D2, state_desc_pint)
         ! ierr = pio_def_var(file, 'state_pintdry',   pio_double, dimids_3D2, state_desc_pintdry)
         ! ierr = pio_def_var(file, 'state_lnpint',    pio_double, dimids_3D2, state_desc_lnpint)
         ! ierr = pio_def_var(file, 'state_lnpintdry', pio_double, dimids_3D2, state_desc_lnpintdry)
         ! ierr = pio_def_var(file, 'state_zi',        pio_double, dimids_3D2, state_desc_zi)
      end if

      !-------------------------------------------------------------------------
      ! define physics tendency variables
      if (add_phys_tend) then
         ! ierr = pio_def_var(file, 'tend_dtdt',   pio_double, dimids_3D1, tend_desc_dtdt )
         ! ierr = pio_def_var(file, 'tend_dudt',   pio_double, dimids_3D1, tend_desc_dudt )
         ! ierr = pio_def_var(file, 'tend_dvdt',   pio_double, dimids_3D1, tend_desc_dvdt)
         ! ierr = pio_def_var(file, 'tend_flx_net',pio_double, dimids_hrz, tend_desc_flx_net)
         ! ierr = pio_def_var(file, 'tend_te_tnd', pio_double, dimids_hrz, tend_desc_te_tnd) ! cumulative boundary flux of total energy
         ! ierr = pio_def_var(file, 'tend_tw_tnd', pio_double, dimids_hrz, tend_desc_tw_tnd) ! cumulative boundary flux of total water
      end if

      !-------------------------------------------------------------------------
      ! define cam_out variables
      if (add_cam_out) then
         ! do m=1,pcnst
         !    write(num,'(i4.4)') m
         !    ierr = pio_def_var(file, 'cam_out_QBOT'//num, pio_double, dimids_hrz, desc_qbot(m))
         ! end do
         ! ierr = pio_def_var(file, 'cam_out_TBOT', pio_double, dimids_hrz, desc_tbot )
         ! ierr = pio_def_var(file, 'cam_out_ZBOT', pio_double, dimids_hrz, desc_zbot )
         ! ierr = pio_def_var(file, 'cam_out_UBOT', pio_double, dimids_hrz, desc_ubot )
         ! ierr = pio_def_var(file, 'cam_out_VBOT', pio_double, dimids_hrz, desc_vbot )
         ! ierr = pio_def_var(file, 'cam_out_PBOT', pio_double, dimids_hrz, desc_pbot )
         ! ierr = pio_def_var(file, 'cam_out_RHO', pio_double, dimids_hrz, desc_rho )
         ierr = pio_def_var(file, 'cam_out_NETSW', pio_double, dimids_hrz, desc_netsw )
         ierr = pio_def_var(file, 'cam_out_FLWDS', pio_double, dimids_hrz, desc_flwds )
         ierr = pio_def_var(file, 'cam_out_PRECSC', pio_double, dimids_hrz, desc_precsc )
         ! ierr = pio_def_var(file, 'cam_out_PRECSL', pio_double, dimids_hrz, desc_precsl )
         ierr = pio_def_var(file, 'cam_out_PRECC', pio_double, dimids_hrz, desc_precc )
         ! ierr = pio_def_var(file, 'cam_out_PRECL', pio_double, dimids_hrz, desc_precl )
         ierr = pio_def_var(file, 'cam_out_SOLS',  pio_double, dimids_hrz, desc_sols )
         ierr = pio_def_var(file, 'cam_out_SOLL',  pio_double, dimids_hrz, desc_soll )
         ierr = pio_def_var(file, 'cam_out_SOLSD', pio_double, dimids_hrz, desc_solsd )
         ierr = pio_def_var(file, 'cam_out_SOLLD', pio_double, dimids_hrz, desc_solld )
         ! ierr = pio_def_var(file, 'cam_out_THBOT', pio_double, dimids_hrz, desc_thbot )
         ! ierr = pio_def_var(file, 'cam_out_PSL', pio_double, dimids_hrz, desc_psl )
         ! ierr = pio_def_var(file, 'cam_out_WSRESP', pio_double, dimids_hrz, desc_wsresp )
         ! ierr = pio_def_var(file, 'cam_out_TAU_EST', pio_double, dimids_hrz, desc_tau_est )
         ! ierr = pio_def_var(file, 'cam_out_UGUST', pio_double, dimids_hrz, desc_ugust )
         ! ierr = pio_def_var(file, 'cam_out_UOVERN', pio_double, dimids_hrz, desc_uovern )
         ! ierr = pio_def_var(file, 'cam_out_CO2PROG', pio_double, dimids_hrz, desc_co2prog )
         ! ierr = pio_def_var(file, 'cam_out_CO2DIAG', pio_double, dimids_hrz, desc_co2diag )
         ! ierr = pio_def_var(file, 'cam_out_BCPHIDRY', pio_double, dimids_hrz, desc_bcphidry )
         ! ierr = pio_def_var(file, 'cam_out_BCPHODRY', pio_double, dimids_hrz, desc_bcphodry )
         ! ierr = pio_def_var(file, 'cam_out_OCPHIDRY', pio_double, dimids_hrz, desc_ocphidry )
         ! ierr = pio_def_var(file, 'cam_out_OCPHODRY', pio_double, dimids_hrz, desc_ocphodry )
         ! ierr = pio_def_var(file, 'cam_out_dstdry1', pio_double, dimids_hrz, desc_dstdry1 )
         ! ierr = pio_def_var(file, 'cam_out_dstdry2', pio_double, dimids_hrz, desc_dstdry2 )
         ! ierr = pio_def_var(file, 'cam_out_dstdry3', pio_double, dimids_hrz, desc_dstdry3 )
         ! ierr = pio_def_var(file, 'cam_out_dstdry4', pio_double, dimids_hrz, desc_dstdry4 )
         ! ierr = pio_def_var(file, 'cam_out_BCPHIWET', pio_double, dimids_hrz, desc_bcphiwet )
         ! ierr = pio_def_var(file, 'cam_out_OCPHIWET', pio_double, dimids_hrz, desc_ocphiwet )
         ! ierr = pio_def_var(file, 'cam_out_dstwet1', pio_double, dimids_hrz, desc_dstwet1 )
         ! ierr = pio_def_var(file, 'cam_out_dstwet2', pio_double, dimids_hrz, desc_dstwet2 )
         ! ierr = pio_def_var(file, 'cam_out_dstwet3', pio_double, dimids_hrz, desc_dstwet3 )
         ! ierr = pio_def_var(file, 'cam_out_dstwet4', pio_double, dimids_hrz, desc_dstwet4 )
      end if

      !-------------------------------------------------------------------------
      ! define cam_in variables
      if (add_cam_in) then
         ierr = pio_def_var(file, 'cam_in_ASDIR',  pio_double, dimids_hrz, desc_asdir)
         ierr = pio_def_var(file, 'cam_in_ASDIF',  pio_double, dimids_hrz, desc_asdif)
         ierr = pio_def_var(file, 'cam_in_ALDIR',  pio_double, dimids_hrz, desc_aldir)
         ierr = pio_def_var(file, 'cam_in_ALDIF',  pio_double, dimids_hrz, desc_aldif)
         ierr = pio_def_var(file, 'cam_in_LWUP',  pio_double, dimids_hrz, desc_lwup)
         ! ierr = pio_def_var(file, 'cam_in_LHF',  pio_double, dimids_hrz, desc_lhf)
         ! ierr = pio_def_var(file, 'cam_in_SHF',  pio_double, dimids_hrz, desc_shf)
         ! ierr = pio_def_var(file, 'cam_in_H2OTEMP',  pio_double, dimids_hrz, desc_h2otemp)
         ! ierr = pio_def_var(file, 'cam_in_WSX',  pio_double, dimids_hrz, desc_wsx)
         ! ierr = pio_def_var(file, 'cam_in_WSY',  pio_double, dimids_hrz, desc_wsy)
         ! ierr = pio_def_var(file, 'cam_in_TREF',  pio_double, dimids_hrz, desc_tref)
         ! ierr = pio_def_var(file, 'cam_in_QREF',  pio_double, dimids_hrz, desc_qref)
         ! ierr = pio_def_var(file, 'cam_in_U10',  pio_double, dimids_hrz, desc_u10)
         ! ierr = pio_def_var(file, 'cam_in_TS',  pio_double, dimids_hrz, desc_ts)
         ! ierr = pio_def_var(file, 'cam_in_SST',  pio_double, dimids_hrz, desc_sst)
         ierr = pio_def_var(file, 'cam_in_SNOWHLAND',  pio_double, dimids_hrz, desc_snowhland)
         ierr = pio_def_var(file, 'cam_in_SNOWHICE',  pio_double, dimids_hrz, desc_snowhice)
         ! ierr = pio_def_var(file, 'cam_in_FCO2_LND',  pio_double, dimids_hrz, desc_fco2_lnd)
         ! ierr = pio_def_var(file, 'cam_in_FCO2_OCN',  pio_double, dimids_hrz, desc_fco2_ocn)
         ! ierr = pio_def_var(file, 'cam_in_FDMS',  pio_double, dimids_hrz, desc_fdms)
         ierr = pio_def_var(file, 'cam_in_LANDFRAC',  pio_double, dimids_hrz, desc_landfrac)
         ierr = pio_def_var(file, 'cam_in_ICEFRAC',  pio_double, dimids_hrz, desc_icefrac)
         ierr = pio_def_var(file, 'cam_in_OCNFRAC',  pio_double, dimids_hrz, desc_ocnfrac)
         ! ierr = pio_def_var(file, 'cam_in_RAM1',  pio_double, dimids_hrz, desc_ram1)
         ! ierr = pio_def_var(file, 'cam_in_FV',  pio_double, dimids_hrz, desc_fv)
         ! ierr = pio_def_var(file, 'cam_in_SOILW',  pio_double, dimids_hrz, desc_soilw) ! MOZART
         ! ierr = pio_def_var(file, 'cam_in_USTAR',  pio_double, dimids_hrz, desc_ustar)
         ! ierr = pio_def_var(file, 'cam_in_RE',  pio_double, dimids_hrz, desc_re)
         ! ierr = pio_def_var(file, 'cam_in_SSQ',  pio_double, dimids_hrz, desc_ssq)

         ! do m=1,pcnst
         !    write(num,'(i4.4)') m
         !    ierr = pio_def_var(file, 'cam_in_CFLX'//num,  pio_double, dimids_hrz, desc_cflx(m))
         ! end do

         ! do m=1,n_drydep
         !    write(num,'(i4.4)') m
         !    ierr = pio_def_var(file, 'cam_in_DEPVEL'//num,  pio_double, dimids_hrz, desc_depvel(m))
         ! end do

         ! do m=1,4
         !    write(num,'(i4.4)') m
         !    ierr = pio_def_var(file, 'cam_in_DSTFLX'//num,  pio_double, dimids_hrz, desc_dstflx(m))
         ! end do

         ! do m=1,shr_megan_mechcomps_n
         !    write(num,'(i4.4)') m
         !    ierr = pio_def_var(file, 'cam_in_MEGANFLX'//num,  pio_double, dimids_hrz, desc_meganflx(m))
         ! end do
      end if

      !-------------------------------------------------------------------------
      ! End variable definitions
      ierr = pio_enddef(file)

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! write data to file
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      
      ! Set missing portions of tmp variables
      ! (only do this once since all chunked vars have the same shape)
      do i=begchunk,endchunk
         if (ncol(i) < pcols) then
            tmp2D(ncol(i)+1:,i)    = fillvalue
            tmp3D(ncol(i)+1:,:,i)  = fillvalue
            tmp3Dp(ncol(i)+1:,:,i) = fillvalue
         end if
      end do

      call cam_grid_write_var(file, grid_id)

      call timemgr_write_restart(file)

      if (add_pbuf) call pbuf_write_restart_alt(file, pbuf2d)

      ! data for prognostic chemistry (does nothing for "chem none" option)
      call chem_write_restart(file)

      !-------------------------------------------------------------------------
      ! write physics state variables
      if (add_phys_state) then
         if (mode==1) then
           do i=begchunk,endchunk
              tmp2D(:ncol(i), i) = phys_state(i)%ps(:ncol(i))
           end do
           call pio_write_darray(file, state_desc_ps, iodesc2d, tmp2D, ierr)
         end if

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = phys_state(i)%psdry(:ncol(i))
         ! end do
         ! call pio_write_darray(file, state_desc_psdry, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = phys_state(i)%phis(:ncol(i))
         ! end do
         ! call pio_write_darray(file, state_desc_phis, iodesc2d, tmp2D, ierr)

         do m=1,3 !m=1,pcnst
            do i=begchunk,endchunk
               tmp3D(:ncol(i),:,i) = phys_state(i)%q(:ncol(i),:,m) 
            end do
            call pio_write_darray(file, state_desc_q(m), iodesc3d, tmp3D, ierr)
         end do

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%t(:ncol(i),:)
         end do
         call pio_write_darray(file, state_desc_t, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%u(:ncol(i),:)
         end do
         call pio_write_darray(file, state_desc_u, iodesc3d, tmp3D, ierr)

         do i=begchunk,endchunk
            tmp3D(:ncol(i),:,i) = phys_state(i)%v(:ncol(i),:) 
         end do
         call pio_write_darray(file, state_desc_v, iodesc3d, tmp3D, ierr)

         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%s(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_s, iodesc3d, tmp3D, ierr)

         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%omega(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_omega, iodesc3d, tmp3D, ierr)
         
         ! state%pmid is only for ml input 
         if (mode==1) then
           do i=begchunk,endchunk
              tmp3D(:ncol(i),:,i) = phys_state(i)%pmid(:ncol(i),:) 
           end do
           call pio_write_darray(file, state_desc_pmid, iodesc3d, tmp3D, ierr)
         end if

         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%pmiddry(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_pmiddry, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%pdel(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_pdel, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%pdeldry(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_pdeldry, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%rpdel(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_rpdel, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%rpdeldry(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_rpdeldry, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%lnpmid(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_lnpmid, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%lnpmiddry(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_lnpmiddry, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%exner(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_exner, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_state(i)%zm(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_zm, iodesc3d, tmp3D, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3Dp(:ncol(i),:,i) = phys_state(i)%pint(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_pint, iodesc3dp, tmp3Dp, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3Dp(:ncol(i),:,i) = phys_state(i)%pintdry(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_pintdry, iodesc3dp, tmp3Dp, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3Dp(:ncol(i),:,i) = phys_state(i)%lnpint(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_lnpint, iodesc3dp, tmp3Dp, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3Dp(:ncol(i),:,i) = phys_state(i)%lnpintdry(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_lnpintdry, iodesc3dp, tmp3Dp, ierr)
         
         ! do i=begchunk,endchunk
         !    tmp3Dp(:ncol(i),:,i) = phys_state(i)%zi(:ncol(i),:) 
         ! end do
         ! call pio_write_darray(file, state_desc_zi, iodesc3dp, tmp3Dp, ierr)

      end if

      !-------------------------------------------------------------------------
      ! write physics state variables
      if (add_phys_tend) then

         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_tend(i)%dtdt(:ncol(i),:)
         ! end do
         ! call pio_write_darray(file, tend_desc_dtdt, iodesc3d, tmp3D, ierr)

         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_tend(i)%dudt(:ncol(i),:)
         ! end do
         ! call pio_write_darray(file, tend_desc_dudt, iodesc3d, tmp3D, ierr)

         ! do i=begchunk,endchunk
         !    tmp3D(:ncol(i),:,i) = phys_tend(i)%dvdt(:ncol(i),:)
         ! end do 
         ! call pio_write_darray(file, tend_desc_dvdt, iodesc3d, tmp3D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i),i) = phys_tend(i)%flx_net(:ncol(i))
         ! end do
         ! call pio_write_darray(file, tend_desc_flx_net, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i),i) = phys_tend(i)%te_tnd(:ncol(i))
         ! end do
         ! call pio_write_darray(file, tend_desc_te_tnd, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i),i) = phys_tend(i)%tw_tnd(:ncol(i))
         ! end do
         ! call pio_write_darray(file, tend_desc_tw_tnd, iodesc2d, tmp2D, ierr)

      end if

      !-------------------------------------------------------------------------
      ! Write cam_in components
      if (add_cam_in) then

         ! do m=1,pcnst
         !    do i=begchunk,endchunk
         !       tmp2D(:ncol(i), i) = cam_in(i)%cflx(:ncol(i), m)
         !    end do
         !    call pio_write_darray(file, desc_cflx(m), iodesc2d, tmp2D, ierr)
         ! end do

         ! do m=1,n_drydep
         !    do i=begchunk,endchunk
         !       tmp2D(:ncol(i), i) = cam_in(i)%depvel(:ncol(i), m)
         !    end do
         !    call pio_write_darray(file, desc_depvel(m), iodesc2d, tmp2D, ierr)
         ! end do

         ! do m=1,4
         !    do i=begchunk,endchunk
         !       tmp2D(:ncol(i), i) = cam_in(i)%dstflx(:ncol(i), m)
         !    end do
         !    call pio_write_darray(file, desc_dstflx(m), iodesc2d, tmp2D, ierr)
         ! end do

         ! do m=1,shr_megan_mechcomps_n
         !    do i=begchunk,endchunk
         !       tmp2D(:ncol(i), i) = cam_in(i)%meganflx(:ncol(i), m)
         !    end do
         !    call pio_write_darray(file, desc_meganflx(m), iodesc2d, tmp2D, ierr)
         ! end do

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%asdir(:ncol(i))
         end do
         call pio_write_darray(file, desc_asdir, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%asdif(:ncol(i))
         end do
         call pio_write_darray(file, desc_asdif, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%aldir(:ncol(i))
         end do
         call pio_write_darray(file, desc_aldir, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%aldif(:ncol(i))
         end do
         call pio_write_darray(file, desc_aldif, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%lwup(:ncol(i))
         end do
         call pio_write_darray(file, desc_lwup, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%lhf(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_lhf, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%shf(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_shf, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%h2otemp(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_h2otemp, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%wsx(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_wsx, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%wsy(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_wsy, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%tref(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_tref, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%qref(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_qref, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%u10(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_u10, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%ts(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_ts, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%sst(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_sst, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%snowhland(:ncol(i))
         end do
         call pio_write_darray(file, desc_snowhland, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%snowhice(:ncol(i))
         end do
         call pio_write_darray(file, desc_snowhice, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%fco2_lnd(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_fco2_lnd, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%fco2_ocn(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_fco2_ocn, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%fdms(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_fdms, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%landfrac(:ncol(i))
         end do
         call pio_write_darray(file, desc_landfrac, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%icefrac(:ncol(i))
         end do
         call pio_write_darray(file, desc_icefrac, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_in(i)%ocnfrac(:ncol(i))
         end do
         call pio_write_darray(file, desc_ocnfrac, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%ram1(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_ram1, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%fv(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_fv, iodesc2d, tmp2D, ierr)

         ! MOZART
         ! do i=begchunk,endchunk
         !   tmp2D(:ncol(i), i) = cam_in(i)%soilw(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_soilw, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%ustar(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_ustar, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%re(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_re, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_in(i)%ssq(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_ssq, iodesc2d, tmp2D, ierr)


      end if

      !-------------------------------------------------------------------------
      ! Write cam_out components
      ! :: tbot(:)     ! bot level temperature
      ! :: zbot(:)     ! bot level height above surface
      ! :: ubot(:)     ! bot level u wind
      ! :: vbot(:)     ! bot level v wind
      ! :: qbot(:,:)   ! bot level specific humidity
      ! :: pbot(:)     ! bot level pressure
      ! :: rho(:)      ! bot level density
      ! :: netsw(:)    !
      ! :: flwds(:)    !
      ! :: precsc(:)   !
      ! :: precsl(:)   !
      ! :: precc(:)    !
      ! :: precl(:)    !
      ! :: soll(:)     !
      ! :: sols(:)     !
      ! :: solld(:)    !
      ! :: solsd(:)    !
      ! :: thbot(:)    !
      ! :: co2prog(:)  ! prognostic co2
      ! :: co2diag(:)  ! diagnostic co2
      ! :: psl(:)
      ! :: wsresp(:)   ! first-order response of low-level wind to surface fluxes
      ! :: tau_est(:)  ! stress estimated to be in equilibrium with ubot/vbot
      ! :: ugust(:)    ! gustiness value
      ! :: uovern(:)       ! ratio of wind speed/brunt vaisalla frequency
            ! (Diagnostic, grid, carbon-flux variables are omitted:
      !  pint, pintdry, lnpint, lnpintdry, zi,
      !  te_ini, te_cur, tw_ini, tw_cur, tc_curr, tc_init, tc_mnst, tc_prev,
      !  c_flux_sfc, c_mflx_air, c_mflx_sff, c_mflx_lnd, c_mflx_ocn,
      !  c_iflx_sfc, c_iflx_air, c_iflx_sff, c_iflx_lnd, c_iflx_ocn
      ! )
      ! :: bcphiwet(:) ! wet deposition of hydrophilic black carbon
      ! :: bcphidry(:) ! dry deposition of hydrophilic black carbon
      ! :: bcphodry(:) ! dry deposition of hydrophobic black carbon
      ! :: ocphiwet(:) ! wet deposition of hydrophilic organic carbon
      ! :: ocphidry(:) ! dry deposition of hydrophilic organic carbon
      ! :: ocphodry(:) ! dry deposition of hydrophobic organic carbon
      ! :: dstwet1(:)  ! wet deposition of dust (bin1)
      ! :: dstdry1(:)  ! dry deposition of dust (bin1)
      ! :: dstwet2(:)  ! wet deposition of dust (bin2)
      ! :: dstdry2(:)  ! dry deposition of dust (bin2)
      ! :: dstwet3(:)  ! wet deposition of dust (bin3)
      ! :: dstdry3(:)  ! dry deposition of dust (bin3)
      ! :: dstwet4(:)  ! wet deposition of dust (bin4)
      ! :: dstdry4(:)  ! dry deposition of dust (bin4)
      if (add_cam_out) then

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%tbot(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_tbot, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%zbot(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_zbot, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%ubot(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_ubot, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%vbot(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_vbot, iodesc2d, tmp2D, ierr)

         ! do m=1,pcnst
         !    do i=begchunk,endchunk
         !       tmp2D(:ncol(i),i) = cam_out(i)%qbot(:ncol(i),m)
         !    end do
         !    call pio_write_darray(file, desc_qbot(m), iodesc2d, tmp2D, ierr)
         ! end do

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%pbot(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_pbot, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%rho(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_rho, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%netsw(:ncol(i))
         end do
         call pio_write_darray(file, desc_netsw, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%flwds(:ncol(i))
         end do
         call pio_write_darray(file, desc_flwds, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%precsc(:ncol(i))
         end do
         call pio_write_darray(file, desc_precsc, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%precsl(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_precsl, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%precc(:ncol(i))
         end do
         call pio_write_darray(file, desc_precc, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%precl(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_precl, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%sols(:ncol(i))
         end do
         call pio_write_darray(file, desc_sols, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%soll(:ncol(i))
         end do
         call pio_write_darray(file, desc_soll, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%solsd(:ncol(i))
         end do
         call pio_write_darray(file, desc_solsd, iodesc2d, tmp2D, ierr)

         do i=begchunk,endchunk
            tmp2D(:ncol(i), i) = cam_out(i)%solld(:ncol(i))
         end do
         call pio_write_darray(file, desc_solld, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%thbot(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_thbot, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%psl(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_psl, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%wsresp(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_wsresp, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%tau_est(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_tau_est, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%ugust(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_ugust, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%uovern(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_uovern, iodesc2d, tmp2D, ierr)

         ! do i=begchunk,endchunk
         !    tmp2D(:ncol(i), i) = cam_out(i)%co2prog(:ncol(i))
         ! end do
         ! call pio_write_darray(file, desc_co2prog, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%co2diag(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_co2diag, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%bcphidry(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_bcphidry, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%bcphodry(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_bcphodry, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%ocphidry(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_ocphidry, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%ocphodry(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_ocphodry, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstdry1(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstdry1, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstdry2(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstdry2, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstdry3(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstdry3, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstdry4(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstdry4, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%bcphiwet(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_bcphiwet, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%ocphiwet(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_ocphiwet, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstwet1(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstwet1, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstwet2(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstwet2, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstwet3(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstwet3, iodesc2d, tmp2D, ierr)

          ! do i=begchunk,endchunk
          !    tmp2D(:ncol(i), i) = cam_out(i)%dstwet4(:ncol(i))
          ! end do
          ! call pio_write_darray(file, desc_dstwet4, iodesc2d, tmp2D, ierr)

      end if

      !-------------------------------------------------------------------------
      ! close the file
      call pio_closefile(file)
      ! call cam_pio_closefile(file)
      
   end subroutine write_ml_training
   !------------------------------------------------------------------------------------------------
end module ml_training
#endif /* MMF_ML_TRAINING */
