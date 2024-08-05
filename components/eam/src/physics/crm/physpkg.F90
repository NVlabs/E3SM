module physpkg
  !-----------------------------------------------------------------------------
  ! Purpose: Provides the interface to MMF physics package
  ! 
  ! Method: 
  !   Each parameterization should be implemented with this sequence of calls:
  !    1)  call the physics routine to calculate tendencies
  !    2)  call physics_update() to apply tendencies from ptend to the state
  !    3)  call check_energy_chng() to ensure that the energy and water 
  !        changes match the boundary fluxes
  !-----------------------------------------------------------------------------
#ifdef SPMD
  use mpishorthand
#endif
  use perf_mod
  use cam_abortutils,          only: endrun
  use shr_kind_mod,            only: i8 => shr_kind_i8, r8 => shr_kind_r8
  use shr_sys_mod,             only: shr_sys_irtc, shr_sys_flush
  use spmd_utils,              only: masterproc, iam
  use physconst,               only: latvap, latice, rh2o
  use physics_types,           only: physics_state, physics_tend, physics_state_set_grid, &
                                     physics_ptend, physics_type_alloc, physics_state_dealloc, physics_tend_dealloc, &
                                     physics_state_alloc, physics_tend_alloc, physics_state_copy
  use physics_update_mod,      only: physics_update, physics_update_init, hist_vars, nvars_prtrb_hist, get_var
  use phys_grid,               only: get_ncols_p, print_cost_p, update_cost_p
  use phys_gmean,              only: gmean_mass
  use ppgrid,                  only: begchunk, endchunk, pcols, pver, pverp
  use constituents,            only: pcnst, cnst_name, cnst_get_ind, setup_moist_indices
  use camsrfexch,              only: cam_out_t, cam_in_t
  use phys_control,            only: phys_do_flux_avg, phys_getopts
  use iop_data_mod,            only: single_column
  use cam_logfile,             only: iulog
  use check_energy,            only: check_energy_set_print_additional_diagn
  implicit none
  private

  !  Physics buffer indices
  integer :: teout_idx            = 0
  integer :: tini_idx             = 0
  integer :: qini_idx             = 0
  integer :: cldliqini_idx        = 0
  integer :: cldiceini_idx        = 0
  integer :: static_ener_ac_idx   = 0
  integer :: water_vap_ac_idx     = 0
  integer :: mmf_clear_rh_idx     = 0
  integer :: species_class(pcnst) = -1 

  save
  !-----------------------------------------------------------------------------
  ! Public methods
  !-----------------------------------------------------------------------------
  public phys_register  ! register physics methods
  public phys_init      ! Public initialization method
  public phys_run1      ! First phase of the public run method
  public phys_run2      ! Second phase of the public run method
  public phys_final     ! Public finalization method
#ifdef MMF_NN_EMULATOR
  public mmf_nn_emulator_driver ! MMF_NN_EMULATOR NN-emulation driver
#endif
  !-----------------------------------------------------------------------------
  ! Private module data
  !-----------------------------------------------------------------------------
  ! Physics package options
  logical :: state_debug_checks   ! Debug physics_state.
  logical :: clim_modal_aero      ! climate controled by prognostic or prescribed modal aerosols
  logical :: prog_modal_aero      ! Prognostic modal aerosols present
  logical :: micro_do_icesupersat
  logical :: pergro_test_active= .false.
  logical :: pergro_mods = .false.
  logical :: is_cmip6_volc       ! true if cmip6 style volcanic file is read otherwise false

contains
!===================================================================================================
!===================================================================================================

subroutine phys_register
  !----------------------------------------------------------------------- 
  ! Purpose: Register constituents and physics buffer fields.
  !-----------------------------------------------------------------------
  use physics_buffer,           only: pbuf_init_time
  use physics_buffer,           only: pbuf_add_field, dtype_r8
  use spmd_utils,               only: masterproc
  use constituents,             only: pcnst, cnst_add, cnst_chk_dim, cnst_name
  use chemistry,                only: chem_register
  use conv_water,               only: conv_water_register
  use physconst,                only: mwdry, cpair, mwh2o, cpwv
  use tracers,                  only: tracers_register
  use check_energy,             only: check_energy_register
  use ghg_data,                 only: ghg_data_register
  use radiation,                only: radiation_register
  use co2_cycle,                only: co2_register
  use co2_diagnostics,          only: co2_diags_register
  use flux_avg,                 only: flux_avg_register
  use ionosphere,               only: ionos_register
  use prescribed_ozone,         only: prescribed_ozone_register
  use prescribed_volcaero,      only: prescribed_volcaero_register
  use prescribed_aero,          only: prescribed_aero_register
  use prescribed_ghg,           only: prescribed_ghg_register
  use sslt_rebin,               only: sslt_rebin_register
  use aoa_tracers,              only: aoa_tracers_register
  use aircraft_emit,            only: aircraft_emit_register
  use cam_diagnostics,          only: diag_register
  use cospsimulator_intr,       only: cospsimulator_intr_register
  use rad_constituents,         only: rad_cnst_get_info ! Added to query if it is a modal aero sim or not
  use output_aerocom_aie,       only: output_aerocom_aie_register, do_aerocom_ind3
  use crm_physics,              only: crm_physics_register
  use vertical_diffusion,       only: vertical_diffusion_register
  use cloud_diagnostics,        only: cloud_diagnostics_register
  use modal_aero_calcsize,      only: modal_aero_calcsize_reg
  use modal_aero_wateruptake,   only: modal_aero_wateruptake_reg
  !---------------------------------------------------------------------------
  ! Local variables
  !---------------------------------------------------------------------------
  integer :: dummy    ! for unused pbuf and constituent indices
  integer :: nmodes
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  call phys_getopts(do_aerocom_ind3_out         = do_aerocom_ind3,  &
                    state_debug_checks_out      = state_debug_checks, &
                    micro_do_icesupersat_out    = micro_do_icesupersat, &
                    pergro_test_active_out      = pergro_test_active, &
                    pergro_mods_out             = pergro_mods )

  ! Initialize dyn_time_lvls
  call pbuf_init_time()

  ! Register water vapor.
  ! This must be the first call to cnst_add so that water vapor is constituent 1.
  call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, dummy, longname='Specific humidity', readiv=.true., is_convtran1=.true.)

  ! Fields for physics package diagnostics
  call pbuf_add_field('TINI',           'physpkg',dtype_r8,(/pcols,pver/), tini_idx)
  call pbuf_add_field('QINI',           'physpkg',dtype_r8,(/pcols,pver/), qini_idx)
  call pbuf_add_field('CLDLIQINI',      'physpkg',dtype_r8,(/pcols,pver/), cldliqini_idx)
  call pbuf_add_field('CLDICEINI',      'physpkg',dtype_r8,(/pcols,pver/), cldiceini_idx)
  call pbuf_add_field('static_ener_ac', 'global', dtype_r8,(/pcols/), static_ener_ac_idx)
  call pbuf_add_field('water_vap_ac',   'global', dtype_r8,(/pcols/), water_vap_ac_idx)
  call pbuf_add_field('vmag_gust',      'global', dtype_r8,(/pcols/), dummy)
  call pbuf_add_field('FRACIS',         'physpkg',dtype_r8,(/pcols,pver,pcnst/),dummy)
  call pbuf_add_field('MMF_CLEAR_RH',   'physpkg',dtype_r8,(/pcols,pver/),mmf_clear_rh_idx)
#ifdef MMF_ML_TRAINING
  call pbuf_add_field('SOLIN',          'global',dtype_r8,(/pcols/), dummy)
  call pbuf_add_field('COSZRS',         'global',dtype_r8,(/pcols/), dummy)
#endif

  ! check energy package
  call check_energy_register()

  ! register fluxes for saving across time
  if (phys_do_flux_avg()) call flux_avg_register()

  call conv_water_register()

  ! Determine whether its a 'modal' aerosol simulation  or not
  call rad_cnst_get_info(0, nmodes=nmodes)
  clim_modal_aero = (nmodes > 0)

  if (clim_modal_aero) then
    call modal_aero_calcsize_reg()
    call modal_aero_wateruptake_reg()
  endif

  ! register chemical constituents including aerosols
  call chem_register(species_class)

  ! co2 constituents
  call co2_register()
  call co2_diags_register()

  call prescribed_volcaero_register()
  call prescribed_ozone_register()
  call prescribed_aero_register()
  call prescribed_ghg_register()
  call sslt_rebin_register

  ! register various data model gasses with pbuf
  call ghg_data_register()

  call aircraft_emit_register()

  call crm_physics_register()

  ! radiation
  call radiation_register()
  call cloud_diagnostics_register()

  ! COSP
  call cospsimulator_intr_register()

  ! vertical diffusion
  call vertical_diffusion_register()

  if (do_aerocom_ind3) call output_aerocom_aie_register()

  ! Register diagnostics PBUF
  call diag_register()

  ! Register age of air tracers
  call aoa_tracers_register()

  ! Register test tracers
  ! This is the last call to register constituents because
  ! the test tracers fill the remaining available slots up
  ! to constituent number PCNST -- regardless of what PCNST is set to.
  call tracers_register()

  ! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()

  ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

end subroutine phys_register

!===================================================================================================
!===================================================================================================

subroutine phys_inidat( cam_out, pbuf2d )
  use physics_buffer,      only: pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, dyn_time_lvls
  use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
  use cam_grid_support,    only: cam_grid_check, cam_grid_id
  use cam_grid_support,    only: cam_grid_get_dim_names
  use pio,                 only: file_desc_t
  use ncdio_atm,           only: infld
  use short_lived_species, only: initialize_short_lived_species
  use comsrf,              only: sgh, sgh30
  use cam_control_mod,     only: aqua_planet
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  integer :: lchnk, m, n, i, k, ncol
  type(file_desc_t), pointer :: fh_ini, fh_topo
  character(len=8) :: fieldname
  real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
  real(r8), pointer :: tptr(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
  real(r8), pointer :: qpert(:,:)

  character*11 :: subname='phys_inidat' ! subroutine name
  character(len=8) :: dim1name, dim2name
  integer :: tpert_idx, qpert_idx, pblh_idx, vmag_gust_idx
  logical :: found=.false., found2=.false.
  integer :: ierr
  integer :: ixcldice, ixcldliq
  integer :: grid_id  ! grid ID for data mapping
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)

  fh_ini=>initial_file_get_id()

  grid_id = cam_grid_id('physgrid')
  if (.not. cam_grid_check(grid_id)) then
    call endrun(trim(subname)//': Internal error, no "physgrid" grid')
  end if
  call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

  if(aqua_planet) then
    sgh = 0._r8
    sgh30 = 0._r8
    if (masterproc) write(iulog,*) 'AQUA_PLANET simulation, sgh, sgh30 initialized to 0.'
  else    
    if (masterproc) write(iulog,*) 'NOT AN AQUA_PLANET simulation, initialize sgh, sgh30, land m using data from file.'
    fh_topo=>topo_file_get_id()
    call infld('SGH', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, sgh, found, gridname='physgrid')
    if(.not. found) call endrun('ERROR: SGH not found on topo file')

    call infld('SGH30', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, sgh30, found, gridname='physgrid')
    if(.not. found) then
      if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
      if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
      sgh30 = sgh
    end if
  end if

  allocate(tptr(1:pcols,begchunk:endchunk))

  call infld('PBLH', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr(:,:), found, gridname='physgrid')
  if(.not. found) then
    tptr(:,:) = 0._r8
    if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
  end if
  pblh_idx = pbuf_get_index('pblh')

  call pbuf_set_field(pbuf2d, pblh_idx, tptr)

  call infld('TPERT', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr(:,:), found, gridname='physgrid')
  if(.not. found) then
    tptr(:,:) = 0._r8
    if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
  end if
  tpert_idx = pbuf_get_index( 'tpert')
  call pbuf_set_field(pbuf2d, tpert_idx, tptr)


  call infld('vmag_gust', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr(:,:), found, gridname='physgrid')
  if(.not. found) then
    tptr(:,:) = 0._r8
    if (masterproc) write(iulog,*) 'vmag_gust initialized to 0.'
  end if
  vmag_gust_idx = pbuf_get_index( 'vmag_gust')
  call pbuf_set_field(pbuf2d, vmag_gust_idx, tptr)


  fieldname='QPERT'  
  qpert_idx = pbuf_get_index( 'qpert',ierr)
  if (qpert_idx > 0) then
    call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr, found, gridname='physgrid')
    if(.not. found) then
      tptr=0_r8
      if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    allocate(tptr3d_2(pcols,pcnst,begchunk:endchunk))
    tptr3d_2 = 0_r8
    tptr3d_2(:,1,:) = tptr(:,:)

    call pbuf_set_field(pbuf2d, qpert_idx, tptr3d_2)
    deallocate(tptr3d_2)
  end if

  fieldname='CUSH'
  m = pbuf_get_index('cush')
  call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr, found, gridname='physgrid')
  if(.not.found) then
    if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
    tptr=1000._r8
  end if
  do n=1,dyn_time_lvls
    call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
  end do
  deallocate(tptr)

  !
  ! 3-D fields
  !

  allocate(tptr3d(pcols,pver,begchunk:endchunk))

  fieldname='CLOUD'
  m = pbuf_get_index('CLD')
  call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if(found) then
    do n = 1, dyn_time_lvls
      call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
    end do
  else
    call pbuf_set_field(pbuf2d, m, 0._r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if

  fieldname='QCWAT'
  m = pbuf_get_index(fieldname,ierr)
  if (m > 0) then
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(.not. found) then
      call infld('Q',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
      if (found) then
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
      else
        call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
      end if
    end if
    do n = 1, dyn_time_lvls
      call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
    end do
  end if

  fieldname = 'ICCWAT'
  m = pbuf_get_index(fieldname, ierr)
  if (m > 0) then
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(found) then
      do n = 1, dyn_time_lvls
        call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
      end do
    else
      call cnst_get_ind('CLDICE', ixcldice)
      call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
      if(found) then
        do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
        end do
      else
        call pbuf_set_field(pbuf2d, m, 0._r8)
      end if
      if (masterproc) then
        if (found) then
          write(iulog,*) trim(fieldname), ' initialized with CLDICE'
        else
          write(iulog,*) trim(fieldname), ' initialized to 0.0'
        end if
      end if
    end if
  end if

  fieldname = 'LCWAT'
  m = pbuf_get_index(fieldname,ierr)
  if (m > 0) then
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(found) then
      do n = 1, dyn_time_lvls
        call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
      end do
    else
      allocate(tptr3d_2(pcols,pver,begchunk:endchunk))     
      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
      call infld('CLDLIQ',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d_2, found2, gridname='physgrid')
      if(found .and. found2) then
        tptr3d(:,:,:)=tptr3d(:,:,:)+tptr3d_2(:,:,:)
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
      else if (found) then ! Data already loaded in tptr3d
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE only'
      else if (found2) then
        tptr3d(:,:,:)=tptr3d_2(:,:,:)
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDLIQ only'
      end if

      if (found .or. found2) then
        do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
        end do
      else
        call pbuf_set_field(pbuf2d, m, 0._r8)
        if (masterproc)  write(iulog,*) trim(fieldname), ' initialized to 0.0'
      end if
      deallocate(tptr3d_2)
    end if
  end if

  deallocate(tptr3d)
  allocate(tptr3d(pcols,pver,begchunk:endchunk))

  fieldname = 'TCWAT'
  m = pbuf_get_index(fieldname,ierr)
  if (m > 0) then
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(.not.found) then
      call infld('T', fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
      if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
    end if
    do n = 1, dyn_time_lvls
      call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
    end do
  end if

  deallocate(tptr3d)
  allocate(tptr3d(pcols,pverp,begchunk:endchunk))

  fieldname = 'TKE'
  m = pbuf_get_index( 'tke')
  call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if (found) then
    call pbuf_set_field(pbuf2d, m, tptr3d)
  else
    call pbuf_set_field(pbuf2d, m, 0.01_r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
  end if


  fieldname = 'KVM'
  m = pbuf_get_index('kvm')
  call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if (found) then
    call pbuf_set_field(pbuf2d, m, tptr3d)
  else
    call pbuf_set_field(pbuf2d, m, 0._r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if


  fieldname = 'KVH'
  m = pbuf_get_index('kvh')
  call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if (found) then
    call pbuf_set_field(pbuf2d, m, tptr3d)
  else
    call pbuf_set_field(pbuf2d, m, 0._r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if

  deallocate(tptr3d)
  allocate(tptr3d(pcols,pver,begchunk:endchunk))

  fieldname = 'CONCLD'
  m = pbuf_get_index('CONCLD')
  call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if(found) then
    do n = 1, dyn_time_lvls
      call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
    end do
  else
    call pbuf_set_field(pbuf2d, m, 0._r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if

  deallocate (tptr3d)

  call initialize_short_lived_species(fh_ini, pbuf2d)

end subroutine phys_inidat

!===================================================================================================
!===================================================================================================

subroutine phys_init( phys_state, phys_tend, pbuf2d, cam_out )
  !-----------------------------------------------------------------------------
  ! Purpose: Initialization of physics package
  !-----------------------------------------------------------------------------
#ifdef MMF_NN_EMULATOR
  use mmf_nn_emulator,            only: init_neural_net
#endif
  use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
  use physconst,          only: rair, cpair, gravit, stebol, tmelt, &
                                latvap, latice, rh2o, rhoh2o, pstd, zvir, &
                                karman, rhodair, physconst_init 
  use ref_pres,           only: pref_edge, pref_mid
  use cloud_rad_props,    only: cloud_rad_props_init
  use cam_control_mod,    only: nsrest  ! restart flag
  use check_energy,       only: check_energy_init, print_additional_diagn
  use chemistry,          only: chem_init
  use prescribed_ozone,   only: prescribed_ozone_init
  use prescribed_ghg,     only: prescribed_ghg_init
  use prescribed_aero,    only: prescribed_aero_init
  use seasalt_model,      only: init_ocean_data, has_mam_mom
  use aerodep_flx,        only: aerodep_flx_init
  use aircraft_emit,      only: aircraft_emit_init
  use prescribed_volcaero,only: prescribed_volcaero_init
  use co2_cycle,          only: co2_init, co2_transport
  use co2_diagnostics,    only: co2_diags_init
  use cam_diagnostics,    only: diag_init
  use gw_drag,            only: gw_init
  use radheat,            only: radheat_init
  use radiation,          only: radiation_init
  use wv_saturation,      only: wv_sat_init
  use ndrop,              only: ndrop_init
  use conv_water,         only: conv_water_init
  use tracers,            only: tracers_init
  use aoa_tracers,        only: aoa_tracers_init
  use rayleigh_friction,  only: rayleigh_friction_init
  use pbl_utils,          only: pbl_utils_init
  use phys_debug_util,    only: phys_debug_init
  use rad_constituents,   only: rad_cnst_init
  use aer_rad_props,      only: aer_rad_props_init
  use sslt_rebin,         only: sslt_rebin_init
  use tropopause,         only: tropopause_init
  use solar_data,         only: solar_data_init
  use rad_solar_var,      only: rad_solar_var_init
  use output_aerocom_aie, only: output_aerocom_aie_init, do_aerocom_ind3
  use dyn_grid,           only: fv_nphys
  use cam_history,        only: addfld, add_default, horiz_only 
  use crm_physics,        only: crm_physics_init 
  use vertical_diffusion, only: vertical_diffusion_init
  use cloud_diagnostics,  only: cloud_diagnostics_init
  use modal_aero_calcsize,   only: modal_aero_calcsize_init
  use modal_aero_wateruptake,only: modal_aero_wateruptake_init
  use nucleate_ice_cam,      only: nucleate_ice_cam_init
  use hetfrz_classnuc_cam,   only: hetfrz_classnuc_cam_init
  use prescribed_macv2,      only: macv2_rad_props_init

  !-----------------------------------------------------------------------------
  ! Input/output arguments
  !-----------------------------------------------------------------------------
  type(physics_state), pointer       :: phys_state(:)
  type(physics_tend ), pointer       :: phys_tend(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)
  !-----------------------------------------------------------------------------
  ! local variables
  !-----------------------------------------------------------------------------
  integer :: lchnk
  real(r8) :: dp1 = huge(1.0_r8) ! set in namelist, assigned in cloud_fraction.F90
  real(r8) :: bulk_scale         ! prescribed aerosol bulk sulfur scale factor
  real(r8), parameter :: mincld = 0.0001_r8 ! minimum allowed cloud fraction for aerosol nucleation
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

  do lchnk = begchunk, endchunk
     call physics_state_set_grid(lchnk, phys_state(lchnk))
  end do

  !-----------------------------------------------------------------------------
  ! Initialize variables in physconst which are not temporally and/or spatially constant
  !----------------------------------------------------------------------------- 
  call physconst_init()

  ! Initialize debugging a physics column
  call phys_debug_init()

  call pbuf_initialize(pbuf2d)

  !initialize physics update interface routine
  call physics_update_init()

  ! diag_init makes addfld calls for dyn fields that are output from the physics decomp
  call diag_init()

  call setup_moist_indices()

  call check_energy_init(phys_state)

  call tracers_init()

  ! age of air tracers
  call aoa_tracers_init()

  teout_idx = pbuf_get_index( 'TEOUT')

  if (nsrest .eq. 0) call phys_inidat(cam_out, pbuf2d) 
  
  ! wv_saturation is relatively independent of everything else and
  ! low level, so init it early. Must at least do this before radiation.
  call wv_sat_init()

  ! Initialize rad constituents and their properties
  call rad_cnst_init()
  call aer_rad_props_init()
  call cloud_rad_props_init()

  ! solar irradiance data modules
  call solar_data_init()

  ! Prognostic chemistry.
  call chem_init(phys_state,pbuf2d, species_class)

  ! Prescribed tracers
  call prescribed_ozone_init()
  call prescribed_ghg_init()
  call prescribed_aero_init()
  call aerodep_flx_init()
  call aircraft_emit_init(phys_state,pbuf2d)
  !when is_cmip6_volc is true ,cmip6 style volcanic file is read
  !Initialized to .false. here but it gets its values from prescribed_volcaero_init
  is_cmip6_volc = .false. 
  call prescribed_volcaero_init(is_cmip6_volc)

  ! Initialize ocean data
  if (has_mam_mom) call init_ocean_data()

  ! Initialize MACv2-SP aerosols
  call macv2_rad_props_init()

  ! co2 cycle            
  if (co2_transport()) call co2_init()
  call co2_diags_init(phys_state)

  call gw_init()

  call rayleigh_friction_init()

  call pbl_utils_init(gravit, karman, cpair, rair, zvir)
  call vertical_diffusion_init(pbuf2d)

  call tsinti(tmelt, latvap, rair, stebol, latice)

  call radiation_init(phys_state)

  call rad_solar_var_init()

  call cloud_diagnostics_init()

  call radheat_init(pref_mid)

  ! The following calls are needed in place of microp_aero_init()
  if (clim_modal_aero) call ndrop_init()
  call nucleate_ice_cam_init(mincld, bulk_scale)
  call hetfrz_classnuc_cam_init(mincld)

  call conv_water_init

  call crm_physics_init(phys_state, pbuf2d, species_class)

  call sslt_rebin_init()
  call tropopause_init()

  if(do_aerocom_ind3) call output_aerocom_aie_init()

  call phys_getopts(prog_modal_aero_out=prog_modal_aero)

  if (clim_modal_aero) then

     ! If climate calculations are affected by prescribed modal aerosols, the
     ! the initialization routine for the dry mode radius calculation is called
     ! here.  For prognostic MAM the initialization is called from
     ! modal_aero_initialize
     if (.not. prog_modal_aero) then
        call modal_aero_calcsize_init(pbuf2d, species_class)
     endif

     call modal_aero_wateruptake_init(pbuf2d)

  end if
  
 !BSINGH - addfld and adddefault calls for perturb growth testing    
  if(pergro_test_active)call add_fld_default_calls()

  !disable additional diagn for crm
  call check_energy_set_print_additional_diagn(.false.)

#ifdef MMF_NN_EMULATOR
  call init_neural_net()
#endif

end subroutine phys_init

!===================================================================================================
!===================================================================================================

#ifdef MMF_NN_EMULATOR
subroutine mmf_nn_emulator_driver(phys_state, phys_state_aphys1, phys_state_mmf, ztodt, phys_tend, pbuf2d,  cam_in, cam_out, cam_out_mmf)
  !-----------------------------------------------------------------------------
  ! Purpose: mmf_nn_emulator driver
  !-----------------------------------------------------------------------------
  use mmf_nn_emulator,          only: cb_partial_coupling, cb_partial_coupling_vars, cb_spinup_step, cb_do_ramp, cb_ramp_linear_steps, &
                              cb_ramp_option, cb_ramp_factor, cb_ramp_step_0steps, cb_ramp_step_1steps
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk, &
                              pbuf_allocate, pbuf_get_index, pbuf_get_field
  use time_manager,     only: get_nstep, get_step_size, & 
                              is_first_step,  is_first_restart_step, &
                              get_curr_calday
  use cam_diagnostics,  only: diag_allocate, &
                              diag_mmf_nn_emulator_debug
  use radiation,        only: iradsw, use_rad_dt_cosz 
  use radconstants,     only: nswbands, get_ref_solar_band_irrad
  use rad_solar_var,    only: get_variability
  use orbit,            only: zenith
  use shr_orb_mod,      only: shr_orb_decl
  use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
  use cam_control_mod,  only: lambm0, obliqr, eccen, mvelpp
  use tropopause,       only: tropopause_output
  use camsrfexch,       only: cam_export
  use cam_diagnostics,  only: diag_export
  use geopotential,        only: geopotential_t
  use cam_history_support, only: pflds
  use physconst,        only: cpair, zvir, rair, gravit
  use string_utils,    only: to_lower
  use cam_diagnostics,        only: diag_phys_writeout, diag_conv
  use cloud_diagnostics,      only: cloud_diagnostics_calc

  !-----------------------------------------------------------------------------
  ! Interface arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_state), intent(in),    dimension(begchunk:endchunk) :: phys_state_aphys1
  type(physics_state), intent(inout),    dimension(begchunk:endchunk) :: phys_state_mmf
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
  type(physics_buffer_desc), pointer, dimension(:,:)               :: pbuf2d
  type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
  type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
  type(cam_out_t),      intent(out),      dimension(begchunk:endchunk) :: cam_out_mmf

  !-----------------------------------------------------------------------------
  ! Local Variables
  !-----------------------------------------------------------------------------
  integer :: lchnk                             ! chunk identifier
  integer :: ncol                              ! number of columns
  integer :: nstep                             ! current timestep number
  type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

  ! stuff for calling crm_physics_tend
  logical           :: use_ECPP
  type(physics_ptend), dimension(begchunk:endchunk) :: ptend ! indivdual parameterization tendencies

  integer  :: itim_old, cldo_idx, cld_idx   ! pbuf indices
  real(r8), pointer, dimension(:,:) :: cld  ! cloud fraction
  real(r8), pointer, dimension(:,:) :: cldo ! old cloud fraction
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Local Variables (for MMF_NN_EMULATOR)
  !-----------------------------------------------------------------------------
  integer :: c, i, j, k 
  integer, save :: nstep0
  integer       :: nstep_NN, dtime
  logical :: do_mmf_nn_emulator_inference = .false.

  ! - for partial coupling - !
! type(physics_state), dimension(begchunk:endchunk)  :: phys_state_nn
  ! type(physics_state), dimension(begchunk:endchunk)  :: phys_state_mmf_backup
  ! type(physics_tend ), dimension(begchunk:endchunk)  :: phys_tend_nn

  type(physics_state), allocatable, dimension(:)  :: phys_state_nn
  type(physics_state), allocatable, dimension(:)  :: phys_state_mmf_backup
  type(physics_tend ), allocatable, dimension(:)  :: phys_tend_nn
  type(cam_out_t),     dimension(begchunk:endchunk)  :: cam_out_nn
  integer :: ixcldice, ixcldliq
  integer :: prec_dp_idx, snow_dp_idx
  real(r8), dimension(:), pointer              :: prec_dp    , snow_dp
  real(r8), dimension(pcols,begchunk:endchunk) :: prec_dp_nn , snow_dp_nn, &
                                                  prec_dp_mmf, snow_dp_mmf
  logical  :: do_geopotential = .false.
  real(r8) :: zvirv_loc(pcols,pver), rairv_loc(pcols,pver)  
  real(r8) :: ramp_ratio 
  ! - !

  real(r8) :: calday       ! current calendar day
  real(r8) :: clat(pcols)  ! current latitudes(radians)
  real(r8) :: clon(pcols)  ! current longitudes(radians)
  real(r8), dimension(pcols,begchunk:endchunk) :: coszrs  ! Cosine solar zenith angle
  real(r8), dimension(pcols,begchunk:endchunk) :: solin   ! Insolation

  real(r8) :: sfac(1:nswbands)  ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band
  real(r8) :: solar_band_irrad(1:nswbands) ! rrtmg-assumed solar irradiance in each sw band
  real(r8) :: dt_avg = 0.0_r8   ! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
  real(r8) :: delta    ! Solar declination angle  in radians
  real(r8) :: eccf     ! Earth orbit eccentricity factor
  integer :: ierr=0
!-----------------------------------------------------------------------------
  ! phys_run1 opening
  ! - phys_timestep_init advances ghg gases,
  ! - need to advance solar insolation (for NN)
  !-----------------------------------------------------------------------------


  allocate(phys_state_nn(begchunk:endchunk), stat=ierr)
   if (ierr /= 0) then
      ! Handle allocation error
      write(iulog,*) 'Error allocating phys_state_nn error = ',ierr
   end if

   allocate(phys_state_mmf_backup(begchunk:endchunk), stat=ierr)
   if (ierr /= 0) then
      ! Handle allocation error
      write(iulog,*) 'Error allocating phys_state_mmf_backup error = ',ierr
   end if

   allocate(phys_tend_nn(begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
        ! Handle allocation error
        write(iulog,*) 'Error allocating phys_tend_nn error = ',ierr
    end if

    do lchnk=begchunk,endchunk
      call physics_state_alloc(phys_state_nn(lchnk),lchnk,pcols)
      call physics_state_alloc(phys_state_mmf_backup(lchnk),lchnk,pcols)
      ! call physics_tend_alloc(phys_tend_nn(lchnk),lchnk,pcols)
   end do

   do lchnk=begchunk,endchunk
      call physics_tend_alloc(phys_tend_nn(lchnk),phys_state_nn(lchnk)%psetcols)
   end do

  nstep = get_nstep()
  dtime = get_step_size()

  call pbuf_allocate(pbuf2d, 'physpkg')
  call diag_allocate()

  ! Advance time information
  call t_startf ('phys_timestep_init')
  call phys_timestep_init( phys_state, cam_out, pbuf2d)
  call t_stopf ('phys_timestep_init')

  ! Calculate  COSZRS and SOLIN
  call get_ref_solar_band_irrad( solar_band_irrad ) ! this can move to init subroutine
  call get_variability(sfac)                        ! "
  do lchnk=begchunk,endchunk
     ncol = phys_state(lchnk)%ncol
     calday = get_curr_calday(-dtime) ! get current calendar day with a negative offset to match the time in mli and CRM physics
     ! coszrs
     call get_rlat_all_p(lchnk, ncol, clat)
     call get_rlon_all_p(lchnk, ncol, clon)
     if (use_rad_dt_cosz)  then
        dtime  = get_step_size()
        dt_avg = iradsw*dtime
     end if
     call zenith(calday, clat, clon, coszrs(:,lchnk), ncol, dt_avg)
     ! solin
     call shr_orb_decl(calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                       delta   ,eccf      )
     solin(:,lchnk) = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(:,lchnk)
  end do
  ! [TO-DO] Check solin and coszrs from this calculation vs. pbuf_XXX

  prec_dp_idx = pbuf_get_index('PREC_DP', errcode=i) ! Query physics buffer index
  snow_dp_idx = pbuf_get_index('SNOW_DP', errcode=i)

  !-----------------------------------------------------------------------------
  ! phys_run1 main
  !-----------------------------------------------------------------------------

  ! Call init subroutine for neural networks
  ! (loading neural network weights and normalization factors)
  if (is_first_step() .or. is_first_restart_step()) then
      ! call init_neural_net()
     nstep0 = nstep
  end if

  ! Determine if MMF spin-up perioid is over
  ! (currently spin up time is set at 86400 sec ~ 1 day)
  ! [TO-DO] create a namelist variable for mmf spin-up time
  ! nstep_NN = 86400 / get_step_size()
  nstep_NN = cb_spinup_step

  if (nstep-nstep0 .eq. nstep_NN) then
     do_mmf_nn_emulator_inference = .true.
     if (masterproc) then
        write(iulog,*) '---------------------------------------'
        write(iulog,*) '[MMF_NN_EMULATOR] NN coupling starts'
        write(iulog,*) '---------------------------------------'
     end if
  end if

#ifdef MMF_NN_EMULATORDEBUG
  if (masterproc) then
     write (iulog,*) '[MMF_NN_EMULATORDEBUG] nstep - nstep0, nstep_NN, do_mmf_nn_emulator = ', nstep - nstep0, nstep_NN, do_mmf_nn_emulator_inference
  endif
#endif

  !Save original values of subroutine arguments
  if (do_mmf_nn_emulator_inference .and. cb_partial_coupling) then
     do lchnk = begchunk, endchunk
      ! since phys_state_mmf_backup etc is just allocated but have not been initialized (empty), doing this copy won't lead to memory leak
        phys_state_nn(lchnk) = phys_state(lchnk) 
        phys_state_mmf_backup(lchnk) = phys_state_mmf(lchnk)
        phys_tend_nn(lchnk)  = phys_tend(lchnk) 
        cam_out_nn(lchnk)    = cam_out(lchnk) 
     end do
  end if

  ! Run phys_run1 physics
  if (.not. do_mmf_nn_emulator_inference) then  ! MMFspin-up
     call phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
     do lchnk = begchunk, endchunk
        call physics_state_dealloc(phys_state_mmf(lchnk)) ! to prevent memory leak
        call physics_state_copy(phys_state(lchnk), phys_state_mmf(lchnk))
        cam_out_mmf(lchnk)    = cam_out(lchnk)
     end do

  else  ! NN inference
     if (cb_partial_coupling) then ! NN partial coupling

        call phys_run1   (phys_state,    ztodt, phys_tend,    pbuf2d, cam_in, cam_out)
        ! store mmf calculation of prec_dp and snow_dp
        do lchnk = begchunk, endchunk
           phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
           call pbuf_get_field(phys_buffer_chunk, prec_dp_idx, prec_dp)
           call pbuf_get_field(phys_buffer_chunk, snow_dp_idx, snow_dp)
           prec_dp_mmf(:,lchnk) = prec_dp(:) 
           snow_dp_mmf(:,lchnk) = snow_dp(:)
        end do

        do lchnk = begchunk, endchunk
          ! update phys_state_mmf but cannot overwrite its the mmf adv and phy history
          call physics_state_dealloc(phys_state_mmf(lchnk))
          call physics_state_copy(phys_state(lchnk), phys_state_mmf(lchnk))
          phys_state_mmf(lchnk)%t_adv(:,:,:) = phys_state_mmf_backup(lchnk)%t_adv(:,:,:)
          phys_state_mmf(lchnk)%u_adv(:,:,:) = phys_state_mmf_backup(lchnk)%u_adv(:,:,:)
          phys_state_mmf(lchnk)%t_phy(:,:,:) = phys_state_mmf_backup(lchnk)%t_phy(:,:,:)
          phys_state_mmf(lchnk)%u_phy(:,:,:) = phys_state_mmf_backup(lchnk)%u_phy(:,:,:)
          phys_state_mmf(lchnk)%q_adv(:,:,:,:) = phys_state_mmf_backup(lchnk)%q_adv(:,:,:,:)
          phys_state_mmf(lchnk)%q_phy(:,:,:,:) = phys_state_mmf_backup(lchnk)%q_phy(:,:,:,:)
          cam_out_mmf(lchnk) = cam_out(lchnk)
        end do

        call phys_run1_NN(phys_state_nn, phys_state_aphys1, ztodt, phys_tend_nn, pbuf2d, cam_in, cam_out_nn,&
                          solin, coszrs)
        ! store nn calculation of prec_dp and snow_dp
        do lchnk = begchunk, endchunk
           phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
           call pbuf_get_field(phys_buffer_chunk, prec_dp_idx, prec_dp)
           call pbuf_get_field(phys_buffer_chunk, snow_dp_idx, snow_dp)
           prec_dp_nn(:,lchnk) = prec_dp(:)
           snow_dp_nn(:,lchnk) = snow_dp(:)
           prec_dp(:) = prec_dp_mmf(:,lchnk) ! restored to mmf calculation
           snow_dp(:) = snow_dp_mmf(:,lchnk) ! (prep for cb_partial_coupling)
        end do

     else ! NN full coupling
        call phys_run1_NN(phys_state, phys_state_aphys1, ztodt, phys_tend, pbuf2d,  cam_in, cam_out,&
                          solin, coszrs)
        do lchnk = begchunk, endchunk
          ! in fully nn coupling case (no partial coupling), phys_state_mmf is just synced with phys_state but won't be used
          call physics_state_dealloc(phys_state_mmf(lchnk))
          call physics_state_copy(phys_state(lchnk), phys_state_mmf(lchnk))
          cam_out_mmf(lchnk)    = cam_out(lchnk)
        end do
     end if ! (cb_partial_coupling)
  end if ! (.not. do_mmf_nn_emulator_inference)

  ! Partial coupling
  ! NN calculations overide MMF calculations for any variables included in 'cb_partial_coupling_vars'
  ! e.g., [ 'ptend_t','ptend_q0001','ptend_q0002','ptend_q0003', 'ptend_u', 'ptend_v',
  !         'cam_out_NETSW', 'cam_out_FLWDS', 'cam_out_PRECSC', 'cam_out_PRECC',
  !         'cam_out_SOLS', 'cam_out_SOLL', 'cam_out_SOLSD', 'cam_out_SOLLD'           ]
  if (do_mmf_nn_emulator_inference .and. cb_partial_coupling) then
  
    if (cb_do_ramp) then
      write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_ramp_option = ', trim(cb_ramp_option)
  
      select case (to_lower(trim(cb_ramp_option)))
        case('constant')
          ramp_ratio = cb_ramp_factor
        case('linear')
          
          if (nstep-nstep0-nstep_NN .le. cb_ramp_linear_steps) then
            ramp_ratio = cb_ramp_factor * (nstep-nstep0-nstep_NN)*1.0/(cb_ramp_linear_steps*1.0)
          else
            ramp_ratio = cb_ramp_factor
          end if
        case('step')
          
          if (mod(nstep-nstep0-nstep_NN, (cb_ramp_step_0steps + cb_ramp_step_1steps)) .le. cb_ramp_step_1steps) then
            ramp_ratio = cb_ramp_factor
          else
            ramp_ratio = 0.0
          end if
      end select
    else
      ramp_ratio = 1.0
    end if
    
    if (cb_do_ramp) then
      write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_do_ramp is on'
      write (iulog,*) 'MMF_NN_EMULATOR partial coupling: 1 is fully NN, ramp_ratio = ', ramp_ratio
      write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_ramp_option = ', trim(cb_ramp_option)
    else
      write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_do_ramp is off'
      write (iulog,*) 'MMF_NN_EMULATOR partial coupling: 1 is fully NN, ramp_ratio = ', ramp_ratio
    end if

     call cnst_get_ind('CLDICE', ixcldice)
     call cnst_get_ind('CLDLIQ', ixcldliq)
     do c = begchunk, endchunk
        k = 1
        do while (k < pflds  .and. cb_partial_coupling_vars(k) /= ' ')
           if (trim(cb_partial_coupling_vars(k)) == 'ptend_t') then
              phys_state(c)%t(:,:)   = phys_state_nn(c)%t(:,:)*ramp_ratio + phys_state(c)%t(:,:)*(1.0_r8-ramp_ratio)
              phys_tend(c)%dtdt(:,:) = phys_tend_nn(c)%dtdt(:,:)*ramp_ratio + phys_tend(c)%dtdt(:,:)*(1.0_r8-ramp_ratio)
              do_geopotential = .true.
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'ptend_q0001') then
              phys_state(c)%q(:,:,1) = phys_state_nn(c)%q(:,:,1)*ramp_ratio + phys_state(c)%q(:,:,1)*(1.0_r8-ramp_ratio) 
              do_geopotential = .true.
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'ptend_q0002') then
              phys_state(c)%q(:,:,ixcldliq) = phys_state_nn(c)%q(:,:,ixcldliq)*ramp_ratio + phys_state(c)%q(:,:,ixcldliq)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'ptend_q0003') then
              phys_state(c)%q(:,:,ixcldice) = phys_state_nn(c)%q(:,:,ixcldice)*ramp_ratio + phys_state(c)%q(:,:,ixcldice)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'ptend_u') then
              phys_state(c)%u(:,:)   = phys_state_nn(c)%u(:,:)*ramp_ratio + phys_state(c)%u(:,:)*(1.0_r8-ramp_ratio)
              phys_tend(c)%dudt(:,:) = phys_tend_nn(c)%dudt(:,:)*ramp_ratio + phys_tend(c)%dudt(:,:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'ptend_v') then
              phys_state(c)%v(:,:)   = phys_state_nn(c)%v(:,:)*ramp_ratio + phys_state(c)%v(:,:)*(1.0_r8-ramp_ratio)
              phys_tend(c)%dvdt(:,:) = phys_tend_nn(c)%dvdt(:,:)*ramp_ratio + phys_tend(c)%dvdt(:,:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_NETSW') then
              cam_out(c)%netsw(:) = cam_out_nn(c)%netsw(:)*ramp_ratio + cam_out(c)%netsw(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_FLWDS') then
              cam_out(c)%flwds(:) = cam_out_nn(c)%flwds(:)*ramp_ratio + cam_out(c)%flwds(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLS') then
              cam_out(c)%sols(:) = cam_out_nn(c)%sols(:)*ramp_ratio + cam_out(c)%sols(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLL') then
              cam_out(c)%soll(:) = cam_out_nn(c)%soll(:)*ramp_ratio + cam_out(c)%soll(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLSD') then
              cam_out(c)%solsd(:) = cam_out_nn(c)%solsd(:)*ramp_ratio + cam_out(c)%solsd(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLLD') then
              cam_out(c)%solld(:) = cam_out_nn(c)%solld(:)*ramp_ratio + cam_out(c)%solld(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_PRECSC') then
              phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
              call pbuf_get_field(phys_buffer_chunk, snow_dp_idx, snow_dp)
              snow_dp(:) = snow_dp_nn(:,c)*ramp_ratio + snow_dp(:)*(1.0_r8-ramp_ratio) 
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_PRECC') then
              phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
              call pbuf_get_field(phys_buffer_chunk, prec_dp_idx, prec_dp)
              prec_dp(:) = prec_dp_nn(:,c)*ramp_ratio + prec_dp(:)*(1.0_r8-ramp_ratio) 
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else
              call endrun('[MMF_NN_EMULATOR: cb_partial_coupling] Wrong variables are included in cb_partial_coupling_vars: ' // trim(cb_partial_coupling_vars(k)))
           end if
           k = k+1
        end do ! k

        if (do_geopotential) then
            ncol = phys_state(c)%ncol
            zvirv_loc(:,:) = zvir
            rairv_loc(:,:) = rair
            call geopotential_t  ( &
                 phys_state(c)%lnpint, phys_state(c)%lnpmid,   phys_state(c)%pint, phys_state(c)%pmid, phys_state(c)%pdel, phys_state(c)%rpdel, &
                 phys_state(c)%t     , phys_state(c)%q(:,:,1), rairv_loc(:,:),  gravit,     zvirv_loc(:,:), &
                 phys_state(c)%zi    , phys_state(c)%zm      , ncol)
            ! update dry static energy for use in next process
            do j = 1, pver
               phys_state(c)%s(:ncol,j) = phys_state(c)%t(:ncol,j)*cpair &
                                          + gravit*phys_state(c)%zm(:ncol,j) + phys_state(c)%phis(:ncol)
            end do ! j
        end if ! (do_geopotential)

     end do ! c
  end if ! (cb_partial coupling)

  ! copy from the tphysbc2 to here. make sure the outputted history file is consistent with the partial coupling
  do lchnk=begchunk, endchunk
    phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
    call t_startf('bc_history_write')
    call diag_phys_writeout(phys_state(lchnk), cam_out(lchnk)%psl)
    call diag_conv(phys_state(lchnk), ztodt, phys_buffer_chunk)
    call t_stopf('bc_history_write')

    !-----------------------------------------------------------------------------
    ! Write cloud diagnostics on history file
    !-----------------------------------------------------------------------------
    call t_startf('bc_cld_diag_history_write')
    call cloud_diagnostics_calc(phys_state(lchnk), phys_buffer_chunk)
    call t_stopf('bc_cld_diag_history_write')
  end do 
  !-----------------------------------------------------------------------------
  ! phys_run1 closing
  ! - tphysbc2 diagnostic (including cam_export)
  !-----------------------------------------------------------------------------
  do lchnk=begchunk, endchunk
     ! Diagnose the location of the tropopause
     call tropopause_output(phys_state(lchnk))

     ! Save atmospheric fields to force surface models
     phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
     call cam_export(phys_state(lchnk), cam_out(lchnk), phys_buffer_chunk)

     ! Write export state to history file
     call diag_export(cam_out(lchnk))
  end do

  do lchnk=begchunk,endchunk
    call physics_state_dealloc(phys_state_nn(lchnk))
    call physics_state_dealloc(phys_state_mmf_backup(lchnk))
    call physics_tend_dealloc(phys_tend_nn(lchnk))
  end do
  deallocate(phys_state_nn)
  deallocate(phys_state_mmf_backup)
  deallocate(phys_tend_nn)

end subroutine mmf_nn_emulator_driver


subroutine phys_run1_NN(phys_state, phys_state_aphys1, ztodt, phys_tend, pbuf2d,  cam_in, cam_out, &
                        solin, coszrs)
  !-----------------------------------------------------------------------------
  ! Purpose: First part of atmos physics before updating of surface components
  !-----------------------------------------------------------------------------
  use mmf_nn_emulator,         only: neural_net, &
                             cb_partial_coupling, cb_partial_coupling_vars
  use physics_buffer,  only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field
  use time_manager,    only: get_nstep
  use check_energy,    only: check_energy_gmean
  use flux_avg,        only: flux_avg_init
  !-----------------------------------------------------------------------------
  ! Interface arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_state), intent(in),    dimension(begchunk:endchunk) :: phys_state_aphys1
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

  type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
  type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
  type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out

  real(r8), intent(in), dimension(pcols,begchunk:endchunk) :: coszrs  ! Cosine solar zenith angle
  real(r8), intent(in), dimension(pcols,begchunk:endchunk) :: solin   ! Insolation

  !-----------------------------------------------------------------------------
  ! Local Variables
  !-----------------------------------------------------------------------------
  type(physics_state)                              :: state
  type(physics_buffer_desc),pointer, dimension(:)  :: phys_buffer_chunk
  type(physics_ptend)                              :: ptend 
  integer :: lchnk                             ! chunk identifier
  integer :: nstep                             ! current timestep number
  integer :: ncol
  integer :: ixcldice, ixcldliq            ! constituent indices for cloud liquid and ice water.
  real(r8), pointer, dimension(:,:) :: tini
  real(r8), pointer, dimension(:,:) :: qini
  real(r8), pointer, dimension(:,:) :: cldliqini
  real(r8), pointer, dimension(:,:) :: cldiceini
  
  nullify(phys_buffer_chunk)
  nullify(tini)
  nullify(qini)
  nullify(cldliqini)
  nullify(cldiceini)

  !-----------------------------------------------------------------------------
  ! phys_run1 opening
  !-----------------------------------------------------------------------------
  nstep = get_nstep()

    ! Set physics tendencies to 0
  do lchnk=begchunk, endchunk
    ncol = phys_state(lchnk)%ncol
    phys_tend(lchnk)%dtdt(:ncol,:pver)  = 0._r8
    phys_tend(lchnk)%dudt(:ncol,:pver)  = 0._r8
    phys_tend(lchnk)%dvdt(:ncol,:pver)  = 0._r8
 end do

  ! The following initialization depends on the import state (cam_in)
  ! being initialized.  This isn't true when cam_init is called, so need
  ! to postpone this initialization to here.
  if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

  ! Compute total energy of input state and previous output state
  call t_startf ('chk_en_gmean')
  call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
  call t_stopf ('chk_en_gmean')

#ifdef TRACER_CHECK
  call gmean_mass ('before tphysbc DRY', phys_state)
#endif

  ! these initial states will be used in tphysac diagnostics
  do lchnk=begchunk, endchunk
     state = phys_state(lchnk)
     phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
     call pbuf_get_field(phys_buffer_chunk, tini_idx, tini)
     call pbuf_get_field(phys_buffer_chunk, qini_idx, qini)
     call pbuf_get_field(phys_buffer_chunk, cldliqini_idx, cldliqini)
     call pbuf_get_field(phys_buffer_chunk, cldiceini_idx, cldiceini)

     ncol = phys_state(lchnk)%ncol
     tini(:ncol,:pver) = state%t(:ncol,:pver)
     call cnst_get_ind('CLDLIQ', ixcldliq)
     call cnst_get_ind('CLDICE', ixcldice)
     qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
     cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
     cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
  end do

  !-----------------------------------------------------------------------------
  ! Neural network
  !-----------------------------------------------------------------------------
  do lchnk=begchunk, endchunk
     phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
     call neural_net (ptend, phys_state(lchnk), phys_state_aphys1(lchnk), &
                      phys_buffer_chunk, &
                      cam_in(lchnk), cam_out(lchnk), &
                      coszrs(:,lchnk), solin(:,lchnk), &
                      ztodt, lchnk)
     call physics_update (phys_state(lchnk), ptend, ztodt, phys_tend(lchnk))
  end do

  !-----------------------------------------------------------------------------
  ! phys_run1 closing
  !-----------------------------------------------------------------------------
#ifdef TRACER_CHECK
  call gmean_mass ('between DRY', phys_state)
#endif

end subroutine phys_run1_NN
#endif /* MMF_NN_EMULATOR */


subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
  !-----------------------------------------------------------------------------
  ! Purpose: First part of atmos physics before updating of surface components
  !-----------------------------------------------------------------------------
  use time_manager,           only: get_nstep
  use cam_diagnostics,        only: diag_allocate, diag_physvar_ic
  use check_energy,           only: check_energy_gmean
  use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, &
                                    pbuf_allocate, pbuf_old_tim_idx, pbuf_get_index
  use comsrf,                 only: fsns, fsnt, flns, sgh, sgh30, flnt, fsds
  use flux_avg,               only: flux_avg_init
  use check_energy,           only: check_energy_chng
  use crm_physics,            only: ncrms, crm_physics_tend
  use crm_ecpp_output_module, only: crm_ecpp_output_type, crm_ecpp_output_initialize, &
                                    crm_ecpp_output_copy, crm_ecpp_output_finalize
  
  !-----------------------------------------------------------------------------
  ! Interface arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

  type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
  type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
  type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
  !-----------------------------------------------------------------------------
  ! Local Variables
  !-----------------------------------------------------------------------------
  integer :: c                                 ! indices
  integer :: ncol                              ! number of columns
  integer :: nstep                             ! current timestep number
  real(r8):: zero(pcols)                       ! array of zeros
#if (! defined SPMD)
  integer :: mpicom = 0
#endif
  integer(i8) :: beg_count                     ! start time for a chunk
  integer(i8) :: end_count                     ! stop time for a chunk
  integer(i8) :: irtc_rate                     ! irtc clock rate
  real(r8)    :: chunk_cost                    ! measured cost per chunk
  type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

  ! stuff for calling crm_physics_tend
  type(physics_ptend), dimension(begchunk:endchunk) :: ptend ! indivdual parameterization tendencies
  logical           :: use_ECPP
  real(r8), dimension(begchunk:endchunk,pcols) :: mmf_qchk_prec_dp  ! CRM precipitation diagostic (liq+ice)  used for check_energy_chng
  real(r8), dimension(begchunk:endchunk,pcols) :: mmf_qchk_snow_dp  ! CRM precipitation diagostic (ice only) used for check_energy_chng
  real(r8), dimension(begchunk:endchunk,pcols) :: mmf_rad_flux      ! CRM radiative flux diagnostic used for check_energy_chng

  type(crm_ecpp_output_type) :: crm_ecpp_output       ! CRM output data for ECPP calculations
  type(crm_ecpp_output_type) :: crm_ecpp_output_chunk ! CRM output data for ECPP calculations (copy)
  integer  :: ncol_sum
  integer  :: icol_beg(begchunk:endchunk)
  integer  :: icol_end(begchunk:endchunk)

  integer  :: itim_old, cldo_idx, cld_idx   ! pbuf indices  
  real(r8), pointer, dimension(:,:) :: cld  ! cloud fraction
  real(r8), pointer, dimension(:,:) :: cldo ! old cloud fraction
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  zero = 0._r8
  nstep = get_nstep()

  call phys_getopts( use_ECPP_out = use_ECPP )
 
  ! The following initialization depends on the import state (cam_in)
  ! being initialized.  This isn't true when cam_init is called, so need
  ! to postpone this initialization to here.
  if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)
 
  ! Compute total energy of input state and previous output state
  call t_startf ('chk_en_gmean')
  call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
  call t_stopf ('chk_en_gmean')

#ifndef MMF_NN_EMULATOR

  call pbuf_allocate(pbuf2d, 'physpkg')
  call diag_allocate()

  !-----------------------------------------------------------------------------
  ! Advance time information
  !-----------------------------------------------------------------------------

  call t_startf ('phys_timestep_init')
  call phys_timestep_init( phys_state, cam_out, pbuf2d)
  call t_stopf ('phys_timestep_init')
#endif

#ifdef TRACER_CHECK
  call gmean_mass ('before tphysbc DRY', phys_state)
#endif

  !-----------------------------------------------------------------------------
  ! Physics tendency before coupler - Phase 1
  !-----------------------------------------------------------------------------

  call t_barrierf('sync_bc_physics', mpicom)
  call t_startf ('bc_physics')
  call t_startf ('bc_physics1')

!$OMP PARALLEL DO PRIVATE (C, beg_count, phys_buffer_chunk, end_count, chunk_cost)
  do c=begchunk, endchunk

    beg_count = shr_sys_irtc(irtc_rate)

    phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

    ! Output physics terms to IC file
    call t_startf ('diag_physvar_ic')
    call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
    call t_stopf ('diag_physvar_ic')

    call tphysbc1(ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), &
                 phys_state(c), phys_tend(c), phys_buffer_chunk, &
                 fsds(1,c), sgh(1,c), sgh30(1,c), &
                 cam_out(c), cam_in(c) )

    end_count = shr_sys_irtc(irtc_rate)
    chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
    call update_cost_p(c, chunk_cost)

  end do

  call t_stopf ('bc_physics1')

  !-----------------------------------------------------------------------------
  ! CRM physics
  !-----------------------------------------------------------------------------  

  if (use_ECPP) then
    call crm_ecpp_output_initialize(crm_ecpp_output_chunk,pcols,pver)
    call crm_ecpp_output_initialize(crm_ecpp_output,      ncrms,pver)
  end if

  call t_startf('crm_physics_tend')
  call crm_physics_tend(ztodt, phys_state, phys_tend, ptend, pbuf2d, cam_in, cam_out, &
                        species_class, crm_ecpp_output, &
                        mmf_qchk_prec_dp, mmf_qchk_snow_dp, mmf_rad_flux )
  call t_stopf('crm_physics_tend')

  do c=begchunk, endchunk
    call physics_update(phys_state(c), ptend(c), ztodt, phys_tend(c))
    call check_energy_chng(phys_state(c), phys_tend(c), "crm_tend", nstep, ztodt, zero, &
                           mmf_qchk_prec_dp(c,:), mmf_qchk_snow_dp(c,:), mmf_rad_flux(c,:))
  end do

  !-----------------------------------------------------------------------------
  ! Physics tendency before coupler - Phase 2
  !-----------------------------------------------------------------------------

  call t_barrierf('sync_bc_physics', mpicom)
  call t_startf ('bc_physics2')

  ! Determine column start and end indices for crm_ecpp_output
  ncol_sum = 0
  do c=begchunk, endchunk
    icol_beg(c) = ncol_sum + 1
    icol_end(c) = ncol_sum + phys_state(c)%ncol
    ncol_sum = ncol_sum + phys_state(c)%ncol
  end do

!$OMP PARALLEL DO PRIVATE (C, beg_count, phys_buffer_chunk, end_count, chunk_cost)
  do c=begchunk, endchunk

    beg_count = shr_sys_irtc(irtc_rate)

    ! Output physics terms to IC file
    phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

    if (use_ECPP) call crm_ecpp_output_copy( crm_ecpp_output, crm_ecpp_output_chunk, icol_beg(c), icol_end(c))

    call tphysbc2(ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), &
                 phys_state(c), phys_tend(c), phys_buffer_chunk, &
                 fsds(1,c), sgh(1,c), sgh30(1,c), &
                 cam_out(c), cam_in(c), crm_ecpp_output_chunk )

    end_count = shr_sys_irtc(irtc_rate)
    chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
    call update_cost_p(c, chunk_cost)

  end do

  call t_stopf ('bc_physics2')
  call t_stopf ('bc_physics')

  if (use_ECPP) then
    call crm_ecpp_output_finalize(crm_ecpp_output_chunk)
    call crm_ecpp_output_finalize(crm_ecpp_output)
  end if

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

#ifdef TRACER_CHECK
  call gmean_mass ('between DRY', phys_state)
#endif

end subroutine phys_run1

!===================================================================================================
!===================================================================================================

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d, cam_out, cam_in )
  !-----------------------------------------------------------------------------
  ! Purpose: Second part of atmos physics after updating of surface models
  !-----------------------------------------------------------------------------
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
  use mo_lightning,     only: lightning_no_prod
  use cam_diagnostics,  only: diag_deallocate, diag_surf
  use comsrf,           only: trefmxav, trefmnav, sgh, sgh30, fsds 
  use physconst,        only: stebol, latvap
  use time_manager,     only: get_nstep
  use check_energy,     only: ieflx_gmean, check_ieflx_fix 
  use phys_control,     only: ieflx_opt
  !-----------------------------------------------------------------------------
  ! Interface arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
  type(physics_buffer_desc), pointer, dimension(:,:)               :: pbuf2d
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
  type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
  !-----------------------------------------------------------------------------
  ! Local Variables
  !-----------------------------------------------------------------------------
  integer :: c                                 ! chunk index
  integer :: ncol                              ! number of columns
  integer :: nstep                             ! current timestep number
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
#if (! defined SPMD)
  integer :: mpicom = 0
#endif
  integer(i8) :: beg_count, end_count, irtc_rate ! for measuring chunk cost
  real(r8):: chunk_cost
  type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk

  call t_barrierf('sync_ac_physics', mpicom)
  call t_startf ('ac_physics')

  nstep = get_nstep()

  !-----------------------------------------------------------------------------
  ! Lightning production of NO
  !-----------------------------------------------------------------------------
  call t_startf ('lightning_no_prod')
  call lightning_no_prod( phys_state, pbuf2d,  cam_in )
  call t_stopf ('lightning_no_prod')

  !-----------------------------------------------------------------------------
  ! calculate the global mean ieflx 
  !-----------------------------------------------------------------------------
  if(ieflx_opt>0) call ieflx_gmean(phys_state, phys_tend, pbuf2d, cam_in, cam_out, nstep)

!$OMP PARALLEL DO PRIVATE (C, beg_count, NCOL, phys_buffer_chunk, end_count, chunk_cost)
  do c=begchunk,endchunk

    beg_count = shr_sys_irtc(irtc_rate)

    ncol = get_ncols_p(c)
    phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

    ! add the implied internal energy flux to sensible heat flux
    if(ieflx_opt>0) call check_ieflx_fix(c, ncol, nstep, cam_in(c)%shf)

    ! surface diagnostics for history files
    call t_startf('diag_surf')
    call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
    call t_stopf('diag_surf')

    call tphysac(ztodt, cam_in(c),  &
        sgh(1,c), sgh30(1,c), cam_out(c),                              &
        phys_state(c), phys_tend(c), phys_buffer_chunk,&
        fsds(1,c))

    end_count = shr_sys_irtc(irtc_rate)
    chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
    call update_cost_p(c, chunk_cost)

  end do ! Chunk loop

  call t_stopf('ac_physics')

#ifdef TRACER_CHECK
  call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

  call pbuf_deallocate(pbuf2d, 'physpkg')

  call pbuf_update_tim_idx()
  call diag_deallocate()

end subroutine phys_run2

!===================================================================================================
!===================================================================================================

subroutine phys_final( phys_state, phys_tend, pbuf2d )
  use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
  use chemistry,      only : chem_final
  use wv_saturation,  only : wv_sat_final
  use crm_physics,    only : crm_physics_final
  use radiation,      only : radiation_final
  !----------------------------------------------------------------------- 
  ! Purpose: Finalization of physics package
  !-----------------------------------------------------------------------
  ! Input/output arguments
  type(physics_state),       pointer :: phys_state(:)
  type(physics_tend ),       pointer :: phys_tend(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  integer :: lchnk
  
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  if(associated(pbuf2d)) then
     call pbuf_deallocate(pbuf2d,'global')
     deallocate(pbuf2d)
  end if


  do lchnk=begchunk,endchunk
    call physics_state_dealloc(phys_state(lchnk))
    call physics_tend_dealloc(phys_tend(lchnk))
 end do

  deallocate(phys_state)
  deallocate(phys_tend)
  
  call t_startf ('chem_final')
  call chem_final
  call t_stopf ('chem_final')

  call t_startf ('wv_sat_final')
  call wv_sat_final
  call t_stopf ('wv_sat_final')

  call t_startf ('radiation_final')
  call radiation_final()
  call t_stopf ('radiation_final')

  call t_startf ('crm_physics_final')
  call crm_physics_final()
  call t_stopf ('crm_physics_final')

  call t_startf ('print_cost_p')
  call print_cost_p
  call t_stopf ('print_cost_p')

end subroutine phys_final

!===================================================================================================
!===================================================================================================

subroutine tphysac (ztodt, cam_in, sgh, sgh30, cam_out, state, tend, pbuf, fsds )
  !----------------------------------------------------------------------------- 
  ! Purpose: 
  !   Tendency physics after coupling to land, sea, and ice models.
  !   Computes the following:
  !     o Radon surface flux and decay (optional)
  !     o Vertical diffusion and planetary boundary layer
  !     o Multiple gravity wave drag calculations (convective, frontal, orographic)
  !-----------------------------------------------------------------------------
  use physics_buffer,     only: physics_buffer_desc, pbuf_set_field, pbuf_get_index, &
                                pbuf_get_field, pbuf_old_tim_idx
  use chemistry,          only: chem_is_active, chem_timestep_tend, chem_emissions
  use cam_diagnostics,    only: diag_phys_tend_writeout, diag_conv
  use gw_drag,            only: gw_tend
  use vertical_diffusion, only: vertical_diffusion_tend
  use rayleigh_friction,  only: rayleigh_friction_tend
  use constituents,       only: cnst_get_ind
  use physics_types,      only: physics_dme_adjust, set_dry_to_wet, physics_state_check
  use tracers,            only: tracers_timestep_tend
  use aoa_tracers,        only: aoa_tracers_timestep_tend
  use physconst,          only: rhoh2o, latvap,latice, rga
  use aero_model,         only: aero_model_drydep
  use check_energy,       only: check_energy_chng, &
                                check_tracers_data, check_tracers_init, &
                                check_tracers_chng, check_tracers_fini
  use time_manager,       only: get_nstep
  use cam_control_mod,    only: aqua_planet 
  use mo_gas_phase_chemdr,only: map2chm
  use clybry_fam,         only: clybry_fam_set
  use charge_neutrality,  only: charge_fix
  use flux_avg,           only: flux_avg_run
  use phys_control,       only: use_qqflx_fixer
  use co2_cycle,          only: co2_cycle_set_ptend

  implicit none
  !-----------------------------------------------------------------------------
  ! Arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
  real(r8), intent(in) :: fsds(pcols)            ! down solar flux
  real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
  real(r8), intent(in) :: sgh30(pcols)           ! Std. deviation of 30s orography for tms

  type(cam_in_t),      intent(inout) :: cam_in
  type(cam_out_t),     intent(inout) :: cam_out
  type(physics_state), intent(inout) :: state
  type(physics_tend ), intent(inout) :: tend
  type(physics_buffer_desc), pointer :: pbuf(:)
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  real(r8), pointer, dimension(:)   :: static_ener_ac_2d  ! Vertically integrated static energy
  real(r8), pointer, dimension(:)   :: water_vap_ac_2d    ! Vertically integrated water vapor
  real(r8), pointer, dimension(:,:) :: tini               !
  real(r8), pointer, dimension(:,:) :: cld                !
  real(r8), pointer, dimension(:,:) :: qini               !
  real(r8), pointer, dimension(:,:) :: cldliqini          !
  real(r8), pointer, dimension(:,:) :: cldiceini          !
  real(r8), pointer, dimension(:,:) :: dtcore             !
  real(r8), pointer, dimension(:,:) :: ast                ! relative humidity cloud fraction 
  type(check_tracers_data):: tracerint      ! tracer mass integrals and cummulative boundary fluxes
  type(physics_ptend) :: ptend              ! indivdual parameterization tendencies
  integer  :: nstep                         ! current timestep number
  real(r8) :: zero(pcols)                   ! array of zeros
  integer  :: lchnk, ncol                   ! chunk identifier
  integer  :: i,k,m                         ! loop iterators
  integer  :: ixcldice, ixcldliq            ! constituent indices for cloud liq and ice
  integer  :: itim_old                      ! pbuf time index
  real(r8) :: obklen(pcols)                 ! Obukhov length for aerosol dry deposition
  real(r8) :: surfric(pcols)                ! surface friction velocity for aerosol dry deposition
  real(r8) :: fh2o(pcols)                   ! h2o flux to balance source from methane chemistry
  real(r8) :: tmp_q     (pcols,pver)        ! tmp variable
  real(r8) :: tmp_cldliq(pcols,pver)        ! tmp variable
  real(r8) :: tmp_cldice(pcols,pver)        ! tmp variable
  real(r8) :: tmp_t     (pcols,pver)        ! tmp variable
  logical  :: l_tracer_aero, l_vdiff, l_rayleigh, l_gw_drag
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  lchnk = state%lchnk
  ncol  = state%ncol
  nstep = get_nstep()
  
  call phys_getopts( l_tracer_aero_out      = l_tracer_aero      &
                    ,l_vdiff_out            = l_vdiff            &
                    ,l_rayleigh_out         = l_rayleigh         &
                    ,l_gw_drag_out          = l_gw_drag          )

  ! Adjust the surface fluxes to reduce instabilities in near sfc layer
  if (phys_do_flux_avg()) call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)

  ! Validate the physics state.
  if (state_debug_checks) call physics_state_check(state, name="before tphysac")

  ! Associate pointers with physics buffer fields
  itim_old = pbuf_old_tim_idx()

  call pbuf_get_field(pbuf, pbuf_get_index('DTCORE'), dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, tini_idx, tini)
  call pbuf_get_field(pbuf, qini_idx, qini)
  call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
  call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

  call pbuf_get_field(pbuf, pbuf_get_index('CLD'),cld,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))
  call pbuf_get_field(pbuf, pbuf_get_index('AST'),ast,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

  !-----------------------------------------------------------------------------
  ! emissions of aerosols and gas-phase chemistry constituents at surface
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) call chem_emissions( state, cam_in )

  !-----------------------------------------------------------------------------
  ! get zero array for energy checker
  !-----------------------------------------------------------------------------
  zero = 0._r8
  call check_tracers_init(state, tracerint)

  !-----------------------------------------------------------------------------
  ! Check if LHF exceeds the total moisture content of the lowest layer
  !-----------------------------------------------------------------------------
  call qneg4('TPHYSAC ', lchnk, ncol, ztodt, &
              state%q(1,pver,1), state%rpdel(1,pver), &
              cam_in%shf, cam_in%lhf, cam_in%cflx )

  !-----------------------------------------------------------------------------
  ! Source/sink terms for advected tracers.
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then

    call t_startf('adv_tracer_src_snk')

    call tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "tracers_timestep_tend", nstep, ztodt, cam_in%cflx)

    call aoa_tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      

    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "aoa_tracers_timestep_tend", nstep, ztodt, cam_in%cflx)

    ! add tendency from aircraft emissions
    call co2_cycle_set_ptend(state, pbuf, ptend)
    call physics_update(state, ptend, ztodt, tend)

    ! Chemistry calculation
    if (chem_is_active()) then
       call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, pbuf,  fh2o, fsds)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
       call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, cam_in%cflx)
    end if ! chem_is_active

    call t_stopf('adv_tracer_src_snk')

  end if ! l_tracer_aero

  !-----------------------------------------------------------------------------
  ! Vertical diffusion/pbl calculation - (pbl, free atmosphere and molecular)
  !-----------------------------------------------------------------------------

  call t_startf('vertical_diffusion_tend')
  call vertical_diffusion_tend(ztodt, state, cam_in%wsx, cam_in%wsy, cam_in%shf, cam_in%cflx, &
                               surfric, obklen, ptend, ast, cam_in%landfrac, sgh30, pbuf )
  call physics_update(state, ptend, ztodt, tend)
  call t_stopf ('vertical_diffusion_tend')
  
  obklen(:) = 0.
  surfric(:) = 0.

  !-----------------------------------------------------------------------------
  ! Rayleigh friction calculation
  !-----------------------------------------------------------------------------
  if (l_rayleigh) then

    call t_startf('rayleigh_friction')
    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('rayleigh_friction')

    call check_energy_chng(state, tend, "vdiff", nstep, ztodt, &
                           cam_in%cflx(:,1), zero, zero, cam_in%shf)
    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

  end if ! l_rayleigh

  !-----------------------------------------------------------------------------
  !  aerosol dry deposition processes
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then

    call t_startf('aero_drydep')
    call aero_model_drydep( state, pbuf, obklen, surfric, cam_in, ztodt, cam_out, ptend )
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('aero_drydep')

    ! enforce charge neutrality
    call charge_fix( ncol, state%q(:,:,:) )

  end if ! l_tracer_aero

  !-----------------------------------------------------------------------------
  ! Gravity wave drag
  !-----------------------------------------------------------------------------
  if (l_gw_drag) then

    call t_startf('gw_tend')
    call gw_tend(state, sgh, pbuf, ztodt, ptend, cam_in)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('gw_tend')
    
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)

  end if ! l_gw_drag

  !-----------------------------------------------------------------------------
  ! Energy budget checks
  !-----------------------------------------------------------------------------
  call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim_old/),(/pcols,1/))       

  tmp_t(:ncol,:pver) = state%t(:ncol,:pver)

  ! store dse after tphysac in buffer
  do k = 1,pver
    dtcore(:ncol,k) = state%t(:ncol,k)
  end do

  call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist

  ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)
  tmp_q     (:ncol,:pver) = state%q(:ncol,:pver,1)
  tmp_cldliq(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
  tmp_cldice(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
  call physics_dme_adjust(state, tend, qini, ztodt)

  call diag_phys_tend_writeout (state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, tmp_t, qini, cldliqini, cldiceini)

  !-----------------------------------------------------------------------------
  ! Set the chlorine and bromine mass mixing ratios
  !-----------------------------------------------------------------------------
  call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

  !-----------------------------------------------------------------------------
  !Integrate column static energy and water vapor
  !-----------------------------------------------------------------------------
  call pbuf_get_field(pbuf, pbuf_get_index('static_ener_ac'), static_ener_ac_2d )
  call pbuf_get_field(pbuf, pbuf_get_index('water_vap_ac'), water_vap_ac_2d )

  static_ener_ac_2d(:) = 0
  water_vap_ac_2d(:)   = 0
  do i = 1,ncol
    do k = 1,pver
      static_ener_ac_2d(i) = static_ener_ac_2d(i) + state%pdel(i,k)*rga*( state%s(i,k) + state%q(i,k,1)*latvap )
      water_vap_ac_2d(i)   = water_vap_ac_2d(i)   + state%pdel(i,k)*rga*  state%q(i,k,1)
    end do
  end do

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call check_tracers_fini(tracerint)

end subroutine tphysac

!===================================================================================================
!===================================================================================================

subroutine tphysbc1(ztodt, fsns, fsnt, flns, flnt, &
                   state, tend, pbuf, fsds, &
                   sgh, sgh30, cam_out, cam_in )
  !----------------------------------------------------------------------------- 
  ! Purpose: Evaluate physics processes BEFORE coupling to sfc components
  !          Phase 1 - energy fixer and dry adjustment
  !
  ! Pass surface fields for separate surface flux calculations
  ! Dump appropriate fields to history file.
  !-----------------------------------------------------------------------------
  use physics_buffer,         only: physics_buffer_desc, pbuf_get_field
  use physics_buffer,         only: pbuf_get_index, pbuf_old_tim_idx
  use physics_buffer,         only: dyn_time_lvls, pbuf_set_field
  use physics_types,          only: physics_ptend_init, physics_ptend_sum, &
                                    physics_state_check, physics_ptend_scale
  use cam_diagnostics,        only: diag_conv_tend_ini, diag_state_b4_phys_write
  use cam_history,            only: outfld, fieldname_len
  use physconst,              only: cpair, latvap, rga
  use constituents,           only: pcnst, qmin, cnst_get_ind
  use time_manager,           only: get_nstep
  use check_energy,           only: check_energy_chng, check_energy_fix, & 
                                    check_water, check_qflx
  use mo_gas_phase_chemdr,    only: map2chm
  use clybry_fam,             only: clybry_fam_adj
  use output_aerocom_aie,     only: do_aerocom_ind3
  use phys_control,           only: use_qqflx_fixer, use_mass_borrower
  use crm_physics,            only: crm_physics_tend, crm_surface_flux_bypass_tend
  use cloud_diagnostics,      only: cloud_diagnostics_calc

  implicit none
  !-----------------------------------------------------------------------------
  ! Interface Arguments
  !-----------------------------------------------------------------------------
  real(r8),            intent(in   ) :: ztodt         ! 2 delta t (model time increment)
  real(r8),            intent(inout) :: fsns(pcols)   ! Surface solar absorbed flux
  real(r8),            intent(inout) :: fsnt(pcols)   ! Net column abs solar flux at model top
  real(r8),            intent(inout) :: flns(pcols)   ! Srf longwave cooling (up-down) flux
  real(r8),            intent(inout) :: flnt(pcols)   ! Net outgoing lw flux at model top
  real(r8),            intent(inout) :: fsds(pcols)   ! Surface solar down flux
  real(r8),            intent(in   ) :: sgh(pcols)    ! Std. deviation of orography
  real(r8),            intent(in   ) :: sgh30(pcols)  ! Std. deviation of 30 s orography for tms
  type(physics_state), intent(inout) :: state
  type(physics_tend ), intent(inout) :: tend
  type(physics_buffer_desc), pointer :: pbuf(:)
  type(cam_out_t),     intent(inout) :: cam_out
  type(cam_in_t),      intent(in)    :: cam_in
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
  type(physics_state)   :: state_alt        ! alt state for CRM input
  integer  :: nstep                         ! current timestep number
  real(r8) :: net_flx(pcols)
  real(r8) :: rtdt                          ! 1./ztodt
  integer  :: lchnk                         ! chunk identifier
  integer  :: ncol                          ! number of atmospheric columns
  integer  :: ierr
  integer  :: i,k,m,ihist                   ! Longitude, level, constituent indices
  integer  :: ixcldice, ixcldliq            ! constituent indices for cloud liquid and ice water.

  ! physics buffer indices
  integer itim_old, ifld

  ! physics buffer fields for total energy and mass adjustment
  real(r8), pointer, dimension(:  ) :: teout
  real(r8), pointer, dimension(:,:) :: tini
  real(r8), pointer, dimension(:,:) :: qini
  real(r8), pointer, dimension(:,:) :: cldliqini
  real(r8), pointer, dimension(:,:) :: cldiceini
  real(r8), pointer, dimension(:,:) :: dtcore
  real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

  ! energy checking variables
  real(r8) :: zero(pcols)                    ! array of zeros
  real(r8) :: flx_heat(pcols)

  logical   :: lq(pcnst)

  character(len=fieldname_len)   :: varname, vsuffix

  real(r8) :: ftem(pcols,pver)         ! tmp space
  real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
  real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor
  real(r8) :: CIDiff(pcols)            ! Difference in vertically integrated static energy

  logical :: l_bc_energy_fix, l_dry_adj

  call phys_getopts( l_bc_energy_fix_out    = l_bc_energy_fix    &
                    ,l_dry_adj_out          = l_dry_adj          )
  
  !-----------------------------------------------------------------------------
  ! Initialize stuff
  !-----------------------------------------------------------------------------
  call t_startf('bc_init')

  zero = 0._r8
  lchnk = state%lchnk
  ncol  = state%ncol
  rtdt = 1._r8/ztodt
  nstep = get_nstep()

  if (pergro_test_active) then 
    !call outfld calls
    do ihist = 1 , nvars_prtrb_hist
      vsuffix  = trim(adjustl(hist_vars(ihist)))
      varname  = trim(adjustl(vsuffix))//'_topphysbc1' ! form variable name
      call outfld( trim(adjustl(varname)),get_var(state,vsuffix), pcols , lchnk )
    enddo
  endif

  call pbuf_get_field(pbuf, pbuf_get_index('static_ener_ac'), static_ener_ac_2d )
  call pbuf_get_field(pbuf, pbuf_get_index('water_vap_ac'), water_vap_ac_2d )

  ! Integrate and compute the difference
  ! CIDiff = difference of column integrated values
  if( nstep == 0 ) then
    CIDiff(:ncol) = 0.0_r8
    call outfld('DTENDTH', CIDiff, pcols, lchnk )
    call outfld('DTENDTQ', CIDiff, pcols, lchnk )
  else
    ! MSE first
    ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
    do k=2,pver
     ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    CIDiff(:ncol) = (ftem(:ncol,1) - static_ener_ac_2d(:ncol))*rtdt

    call outfld('DTENDTH', CIDiff, pcols, lchnk )
    ! Water vapor second
    ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
    do k=2,pver
      ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    CIDiff(:ncol) = (ftem(:ncol,1) - water_vap_ac_2d(:ncol))*rtdt

    call outfld('DTENDTQ', CIDiff, pcols, lchnk )
  end if

  ! Associate pointers with physics buffer fields
  itim_old = pbuf_old_tim_idx()
  call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))
  call pbuf_get_field(pbuf, tini_idx, tini)
  call pbuf_get_field(pbuf, qini_idx, qini)
  call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
  call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)
  call pbuf_get_field(pbuf, pbuf_get_index('DTCORE'), dtcore, (/1,1,itim_old/), (/pcols,pver,1/) )
  call pbuf_get_field(pbuf, pbuf_get_index('FRACIS'), fracis, (/1,1,1/),        (/pcols, pver, pcnst/)  )
  fracis (:ncol,:,1:pcnst) = 1._r8

  ! Set physics tendencies to 0
  tend%dTdt(:ncol,:pver)  = 0._r8
  tend%dudt(:ncol,:pver)  = 0._r8
  tend%dvdt(:ncol,:pver)  = 0._r8

  !-----------------------------------------------------------------------------
  ! Mass checks and fixers
  !-----------------------------------------------------------------------------

  call check_qflx (state, tend, "PHYBC01", nstep, ztodt, cam_in%cflx(:,1))
  call check_water(state, tend, "PHYBC01", nstep, ztodt)

  ! make sure tracers are all positive - if use_mass_borrower then just print diagnostics
  call qneg3('TPHYSBCb', lchnk, ncol, pcols, pver, 1, pcnst, qmin, state%q, .not.use_mass_borrower )

  if(use_mass_borrower) then 
    ! tracer borrower for mass conservation 
     do m = 1, pcnst 
        call massborrow("PHYBC01",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
     end do
  end if 

  ! Validate state coming from the dynamics.
  if (state_debug_checks) call physics_state_check(state, name="before tphysbc (dycore?)")

  ! Adjust chemistry for conservation issues
  call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

  ! Validate output of clybry_fam_adj
  if (state_debug_checks) call physics_state_check(state, name="clybry_fam_adj")

  ! make sure tracers are all positive, again - if use_mass_borrower then just print diagnostics
  call qneg3('TPHYSBCc',lchnk  ,ncol, pcols, pver, 1, pcnst, qmin, state%q, .not.use_mass_borrower )

  if(use_mass_borrower) then
     ! tracer borrower for mass conservation 
     do m = 1, pcnst
        call massborrow("PHYBC02",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
     end do
  end if

  call check_water(state, tend, "PHYBC02", nstep, ztodt)

  ! Dump out "before physics" state
  call diag_state_b4_phys_write (state)

  call t_stopf('bc_init')

  !-----------------------------------------------------------------------------
  ! Global mean total energy fixer
  !-----------------------------------------------------------------------------
  if (l_bc_energy_fix) then

    call t_startf('energy_fixer')

    tini(:ncol,:pver) = state%t(:ncol,:pver)
    
    call check_energy_fix(state, ptend, nstep, flx_heat)
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
    
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)

    qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
    cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

    call outfld('TEOUT', teout       , pcols, lchnk   )
    call outfld('TEINP', state%te_ini, pcols, lchnk   )
    call outfld('TEFIX', state%te_cur, pcols, lchnk   )

    ! set and output the dse change due to dynpkg
    if( nstep > dyn_time_lvls-1 ) then
      do k = 1,pver
        dtcore(:ncol,k) = (tini(:ncol,k) - dtcore(:ncol,k))/(ztodt) + tend%dTdt(:ncol,k)
      end do
      call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')

  end if
    
  !-----------------------------------------------------------------------------
  ! Dry adjustment
  !-----------------------------------------------------------------------------
  if (l_dry_adj) then

    call t_startf('dry_adjustment')

    ! Copy state info for input to dadadj
    ! This is a kludge so dadadj doesn't have to be reformulated for DSE
    ! This code block is not a good example of interfacing a parameterization
    lq(:) = .FALSE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
    ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
    ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

    call dadadj (lchnk, ncol, state%pmid, state%pint, state%pdel, ptend%s, ptend%q(1,1,1))

    ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
    ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt

    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('dry_adjustment')

  end if

  !-----------------------------------------------------------------------------
  ! MMF surface flux bypass
  !-----------------------------------------------------------------------------
#if defined( MMF_FLUX_BYPASS )

  ! Check if LHF exceeds the total moisture content of the lowest layer
  call qneg4('TPHYSBC ', lchnk, ncol, ztodt, &
              state%q(1,pver,1), state%rpdel(1,pver), &
              cam_in%shf, cam_in%lhf, cam_in%cflx )

  call crm_surface_flux_bypass_tend(state, cam_in, ptend)
  call physics_update(state, ptend, ztodt, tend)  
  call check_energy_chng(state, tend, "crm_tend", nstep, ztodt,  &
                         cam_in%shf(:), zero, zero, cam_in%cflx(:,1)) 
#endif

end subroutine tphysbc1

!===================================================================================================
!===================================================================================================

subroutine tphysbc2(ztodt, fsns, fsnt, flns, flnt, &
                   state, tend, pbuf, fsds, &
                   sgh, sgh30, cam_out, cam_in, crm_ecpp_output )
  !----------------------------------------------------------------------------- 
  ! Purpose: Evaluate physics processes BEFORE coupling to sfc components
  !          Phase 2 - aerosols, radiation, and diagnostics
  !
  ! Pass surface fields for separate surface flux calculations
  ! Dump appropriate fields to history file.
  !-----------------------------------------------------------------------------
  use physics_buffer,         only: physics_buffer_desc, pbuf_get_field
  use physics_buffer,         only: pbuf_get_index, pbuf_old_tim_idx
  use physics_buffer,         only: dyn_time_lvls
  use physics_types,          only: physics_ptend_init, physics_ptend_sum, &
                                    physics_state_check, physics_ptend_scale
  use cam_diagnostics,        only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, &
                                    diag_export, diag_state_b4_phys_write
  use cam_history,            only: outfld, fieldname_len
  use physconst,              only: cpair, latvap, gravit, rga
  use constituents,           only: pcnst, qmin, cnst_get_ind
  use time_manager,           only: get_nstep
  use check_energy,           only: check_energy_chng, & 
                                    check_tracers_data, check_tracers_init, &
                                    check_tracers_chng, check_tracers_fini
  use aero_model,             only: aero_model_wetdep
  use modal_aero_calcsize,    only: modal_aero_calcsize_sub
  use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
  use radiation,              only: radiation_tend
  use tropopause,             only: tropopause_output
  use output_aerocom_aie,     only: do_aerocom_ind3, cloud_top_aerocom
  use cloud_diagnostics,      only: cloud_diagnostics_calc
  use crm_ecpp_output_module, only: crm_ecpp_output_type
  use camsrfexch,             only: cam_export
#if defined( ECPP )
   use module_ecpp_ppdriver2, only: parampollu_driver2
   use module_data_ecpp1,     only: dtstep_pp_input
   use crmclouds_camaerosols, only: crmclouds_mixnuc_tend
#endif
  implicit none
  !-----------------------------------------------------------------------------
  ! Interface Arguments
  !-----------------------------------------------------------------------------
  real(r8),                  intent(in   ) :: ztodt         ! 2 delta t (model time increment)
  real(r8),                  intent(inout) :: fsns(pcols)   ! Surface solar absorbed flux
  real(r8),                  intent(inout) :: fsnt(pcols)   ! Net column abs solar flux at model top
  real(r8),                  intent(inout) :: flns(pcols)   ! Srf longwave cooling (up-down) flux
  real(r8),                  intent(inout) :: flnt(pcols)   ! Net outgoing lw flux at model top
  real(r8),                  intent(inout) :: fsds(pcols)   ! Surface solar down flux
  real(r8),                  intent(in   ) :: sgh(pcols)    ! Std. deviation of orography
  real(r8),                  intent(in   ) :: sgh30(pcols)  ! Std. deviation of 30 s orography for tms
  type(physics_state),       intent(inout) :: state
  type(physics_tend ),       intent(inout) :: tend
  type(physics_buffer_desc), pointer       :: pbuf(:)
  type(cam_out_t),           intent(inout) :: cam_out
  type(cam_in_t),            intent(in   ) :: cam_in
  type(crm_ecpp_output_type),intent(inout) :: crm_ecpp_output   ! CRM output data for ECPP calculations
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
  type(physics_state)   :: state_alt        ! alt state for CRM input
  integer  :: nstep                         ! current timestep number
  real(r8) :: net_flx(pcols)
  real(r8) :: cmfmc2(pcols,pverp)           ! Moist convection cloud mass flux
  real(r8) :: dlf(pcols,pver)               ! Detraining cld H20 from shallow + deep convections
  real(r8) :: dlf2(pcols,pver)              ! Detraining cld H20 from shallow convections
  integer  :: lchnk                         ! chunk identifier
  integer  :: ncol                          ! number of atmospheric columns
  integer  :: ierr
  integer  :: i,k,m,ihist                   ! Longitude, level, constituent indices

  real(r8) :: sh_e_ed_ratio(pcols,pver)      ! shallow conv [ent/(ent+det)] ratio 

  ! pbuf fields
  integer itim_old
  real(r8), pointer, dimension(:,:) :: cld  ! cloud fraction
  real(r8), pointer, dimension(:,:) :: cldo ! old cloud fraction

  ! energy checking variables
  real(r8) :: zero(pcols)                    ! array of zeros
  type(check_tracers_data):: tracerint       ! energy integrals and cummulative boundary fluxes
  real(r8) :: zero_tracers(pcols,pcnst)

  logical   :: lq(pcnst)

  character(len=fieldname_len)   :: varname, vsuffix

  !BSINGH - these were moved from zm_conv_intr because they are  used by aero_model_wetdep 
  real(r8), dimension(pcols,pver) :: mu, eu, du, md, ed, dp
  real(r8):: dsubcld(pcols) ! wg layer thickness in mbs (between upper/lower interface).
  integer :: jt(pcols)      ! wg layer thickness in mbs between lcl and maxi.    
  integer :: maxg(pcols)    ! wg top  level index of deep cumulus convection.
  integer :: ideep(pcols)   ! wg gathered values of maxi.
  integer :: lengath        ! w holds position of gathered points vs longitude index

  logical :: l_tracer_aero, l_rad
  
  logical           :: use_ECPP
  real(r8), pointer :: mmf_clear_rh(:,:) ! CRM clear air relative humidity used for aerosol water uptake

  ! ECPP variables
  real(r8),pointer,dimension(:)   :: pblh              ! PBL height (for ECPP)
  real(r8),pointer,dimension(:,:) :: acldy_cen_tbeg    ! cloud fraction
  real(r8)                        :: dtstep_pp         ! ECPP time step (seconds)
  integer                         :: necpp             ! number of GCM time steps in which ECPP is called once

  call phys_getopts( l_tracer_aero_out      = l_tracer_aero      &
                    ,l_rad_out              = l_rad              &
                    ,use_ECPP_out           = use_ECPP           )
  
  !-----------------------------------------------------------------------------
  ! Initialize stuff
  !-----------------------------------------------------------------------------
  call t_startf('bc_init')

  zero = 0._r8
  zero_tracers(:,:) = 0._r8
  lchnk = state%lchnk
  ncol  = state%ncol
  nstep = get_nstep()

  if (pergro_test_active) then 
    !call outfld calls
    do ihist = 1 , nvars_prtrb_hist
      vsuffix  = trim(adjustl(hist_vars(ihist)))
      varname  = trim(adjustl(vsuffix))//'_topphysbc2' ! form variable name
      call outfld( trim(adjustl(varname)),get_var(state,vsuffix), pcols , lchnk )
    enddo
  endif

  call pbuf_get_field(pbuf, mmf_clear_rh_idx, mmf_clear_rh )

  ! compute mass integrals of input tracers state
  call check_tracers_init(state, tracerint)

  !-----------------------------------------------------------------------------
  ! Modal aerosol wet radius for radiative calculation
  !-----------------------------------------------------------------------------
#if defined( ECPP ) && defined(MODAL_AERO)
  ! temporarily turn on all lq, so it is allocated
  lq(:) = .true.
  call physics_ptend_init(ptend, state%psetcols, 'crm - modal_aero_wateruptake_dr', lq=lq)

  call t_startf('modal_aero_mmf')

  ! set all ptend%lq to false as they will be set in modal_aero_calcsize_sub
  ptend%lq(:) = .false.
  call modal_aero_calcsize_sub (state, ztodt, pbuf, ptend)
  call modal_aero_wateruptake_dr(state, pbuf, clear_rh_in=mmf_clear_rh)

  ! ECPP handles aerosol wet deposition, so tendency from wet depostion is 
  ! not updated in mz_aero_wet_intr (mz_aerosols_intr.F90), but tendencies
  ! from other parts of crmclouds_aerosol_wet_intr() are still updated here.
  call physics_update (state, ptend, ztodt, tend)
  call t_stopf('modal_aero_mmf')

  call check_energy_chng(state, tend, "modal_aero_mmf", nstep, ztodt, &
                         zero, zero, zero, zero)

#endif /* ECPP and MODAL_AERO */
  !-----------------------------------------------------------------------------
  ! ECPP - Explicit-Cloud Parameterized-Pollutant
  !-----------------------------------------------------------------------------
#if defined( ECPP )
  if (use_ECPP) then

    call pbuf_get_field(pbuf, pbuf_get_index('pblh'), pblh)
    call pbuf_get_field(pbuf, pbuf_get_index('ACLDY_CEN'), acldy_cen_tbeg)

    dtstep_pp = dtstep_pp_input
    necpp = dtstep_pp/ztodt

    if (nstep.ne.0 .and. mod(nstep, necpp).eq.0) then

      ! aerosol tendency from droplet activation and mixing
      ! cldo and cldn are set to be the same in crmclouds_mixnuc_tend,
      ! So only turbulence mixing is done here.
      call t_startf('crmclouds_mixnuc')
      call crmclouds_mixnuc_tend(state, ptend, dtstep_pp,           &
                                 cam_in%cflx, pblh, pbuf,           &
                                 crm_ecpp_output%wwqui_cen,         &
                                 crm_ecpp_output%wwqui_cloudy_cen,  &
                                 crm_ecpp_output%wwqui_bnd,         &
                                 crm_ecpp_output%wwqui_cloudy_bnd,  &
                                 species_class)
      call physics_update(state, ptend, dtstep_pp, tend)
      call t_stopf('crmclouds_mixnuc')

      ! ECPP interface
      call t_startf('ecpp')
      call parampollu_driver2(state, ptend, pbuf, dtstep_pp, dtstep_pp,   &
                              crm_ecpp_output%acen,       crm_ecpp_output%abnd,         &
                              crm_ecpp_output%acen_tf,    crm_ecpp_output%abnd_tf,      &
                              crm_ecpp_output%massflxbnd, crm_ecpp_output%rhcen,        &
                              crm_ecpp_output%qcloudcen,  crm_ecpp_output%qlsinkcen,    &
                              crm_ecpp_output%precrcen,   crm_ecpp_output%precsolidcen, &
                              acldy_cen_tbeg)
      call physics_update(state, ptend, dtstep_pp, tend)
      call t_stopf ('ecpp')

    end if ! nstep.ne.0 .and. mod(nstep, necpp).eq.0

  end if ! use_ECPP
#endif /* ECPP */
  !-----------------------------------------------------------------------------
  ! save old CRM cloud fraction - w/o CRM, this is done in cldwat2m.F90
  !-----------------------------------------------------------------------------
  itim_old = pbuf_old_tim_idx()
  call pbuf_get_field(pbuf,pbuf_get_index('CLDO'),cldo,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))
  call pbuf_get_field(pbuf,pbuf_get_index('CLD') ,cld ,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

  cldo(1:ncol,1:pver) = cld(1:ncol,1:pver)

  !-----------------------------------------------------------------------------
  ! Aerosol stuff
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then
    if (use_ECPP) then
      ! With MMF + ECPP we can skip the conventional aerosol routines
    else
      ! Aerosol wet chemistry determines scavenging and transformations.
      ! This is followed by convective transport of all trace species except
      ! water vapor and condensate. Scavenging needs to occur prior to
      ! transport in order to determine interstitial fraction.

      ! Without ECPP we should be using prescribed aerosols, so we only
      ! need to consider the wet deposition and water uptake for radiation

      ! Aerosol wet removal (including aerosol water uptake)
      call t_startf('aero_model_wetdep')
      call aero_model_wetdep( ztodt, dlf, dlf2, cmfmc2, state,  & ! inputs
              sh_e_ed_ratio, mu, md, du, eu, ed, dp, dsubcld,    &
              jt, maxg, ideep, lengath, species_class,           &
              cam_out, pbuf, ptend,                              & ! outputs
              clear_rh=mmf_clear_rh) ! clear air relative humidity for water uptake
      call physics_update(state, ptend, ztodt, tend)
      call t_stopf('aero_model_wetdep')

      ! check tracer integrals
      call check_tracers_chng(state, tracerint, "aero_model_wetdep", nstep, ztodt,  zero_tracers)

    end if
  end if ! l_tracer_aero

#ifndef MMF_NN_EMULATOR
! in mmf_nn_emulator_driver, this block is included, so need to skip when MMF_NN_EMULATOR is defined
! even not executing this block, it won't influence the radiation_tend, which won't use cam_out%psl
  !-----------------------------------------------------------------------------
  ! Moist physical parameteriztions complete: 
  ! send dynamical variables, and derived variables to history file
  !-----------------------------------------------------------------------------
  call t_startf('bc_history_write')
  call diag_phys_writeout(state, cam_out%psl)
  call diag_conv(state, ztodt, pbuf)
  call t_stopf('bc_history_write')

  !-----------------------------------------------------------------------------
  ! Write cloud diagnostics on history file
  !-----------------------------------------------------------------------------
  call t_startf('bc_cld_diag_history_write')
  call cloud_diagnostics_calc(state, pbuf)
  call t_stopf('bc_cld_diag_history_write')
#endif

  !-----------------------------------------------------------------------------
  ! Radiation computations
  !-----------------------------------------------------------------------------
  if (l_rad) then
    call t_startf('radiation')
    call radiation_tend(state,ptend, pbuf, cam_out, cam_in, &
                        cam_in%landfrac, cam_in%icefrac, cam_in%snowhland, &
                        fsns, fsnt, flns, flnt, fsds, &
                        net_flx, is_cmip6_volc, ztodt, clear_rh=mmf_clear_rh)
    call t_stopf('radiation')

    ! We don't need to call physics_update or check_energy_chng for the MMF
    ! because the radiative tendency is added within the call to crm_physics_tend
    ! seems that it will update the crm_rad etc in pbuf, this will update the state in the next call to crm_physics_tend and crm_module

  end if ! l_rad

  !-----------------------------------------------------------------------------
  ! Diagnostics
  !-----------------------------------------------------------------------------

  call t_startf('tphysbc_diagnostics')

  if(do_aerocom_ind3) call cloud_top_aerocom(state, pbuf) 

#ifndef MMF_NN_EMULATOR
! in mmf_nn_emulator_driver, this block is included, so need to skip when MMF_NN_EMULATOR is defined
  ! Diagnose the location of the tropopause
  call tropopause_output(state)

  ! Save atmospheric fields to force surface models
  call cam_export(state,cam_out,pbuf)
  
  ! Write export state to history file
  call diag_export(cam_out)
#endif

  call check_tracers_fini(tracerint)

  call t_stopf('tphysbc_diagnostics')

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
end subroutine tphysbc2

!===================================================================================================
!===================================================================================================

subroutine phys_timestep_init(phys_state, cam_out, pbuf2d)
  !-----------------------------------------------------------------------------
  ! Purpose: The place for parameterizations to call per timestep initializations.
  !          This is used to update time interpolated fields from boundary datasets.
  !-----------------------------------------------------------------------------
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use physics_buffer,      only: physics_buffer_desc
  use ghg_data,            only: ghg_data_timestep_init
  use radiation,           only: radiation_do
  use tracers,             only: tracers_timestep_init
  use aoa_tracers,         only: aoa_tracers_timestep_init
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
  use efield,              only: get_efield
  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv
  use seasalt_model,       only: advance_ocean_data, has_mam_mom

  implicit none
  !-----------------------------------------------------------------------------
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out

  type(physics_buffer_desc), pointer                 :: pbuf2d(:,:)
  !-----------------------------------------------------------------------------

  ! Chemistry surface values
  call chem_surfvals_set()

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call chem_timestep_init(phys_state, pbuf2d)

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state, pbuf2d)
  call prescribed_ghg_adv(phys_state, pbuf2d)
  call prescribed_aero_adv(phys_state, pbuf2d)
  call aircraft_emit_adv(phys_state, pbuf2d)
  call prescribed_volcaero_adv(phys_state, pbuf2d)

  ! advance ocean data for sea salt emissions
  if (has_mam_mom) call advance_ocean_data(phys_state, pbuf2d)

  ! prescribed aerosol deposition fluxes
  call aerodep_flx_adv(phys_state, pbuf2d, cam_out)

  ! Time interpolate data models of gasses in pbuf2d
  call ghg_data_timestep_init(pbuf2d,  phys_state)

  ! Upper atmosphere radiative processes
  call radheat_timestep_init(phys_state, pbuf2d)
 
  ! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init(pbuf2d, phys_state)

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

end subroutine phys_timestep_init

!===================================================================================================
!===================================================================================================

subroutine add_fld_default_calls()
  use cam_history,        only: addfld, add_default, fieldname_len

  implicit none

  !Add all existing ptend names for the addfld calls
  character(len=20), parameter :: vlist(11) = (/'topphysbc           ','chkenergyfix        ',&
                                                'dadadj              ','pcwdetrain_mac      ',&
                                                'aero_model_wetdep_ma','cam_radheat         ',&
                                                'chemistry           ','vdiff               ',&
                                                'rayleigh_friction   ','aero_model_drydep_ma',&
                                                'Grav_wave_drag      ' /)

  character(len=fieldname_len) :: varname
  character(len=1000)          :: s_lngname,stend_lngname,qv_lngname,qvtend_lngname,t_lngname
  
  integer :: iv, ntot, ihist

  ntot = size(vlist)

  do ihist = 1 , nvars_prtrb_hist
     do iv = 1, ntot   
        
        varname  = trim(adjustl(hist_vars(ihist)))//'_'//trim(adjustl(vlist(iv))) ! form variable name

        ! The units and longname are dummy as it is for a test only
        call addfld (trim(adjustl(varname)), (/ 'lev' /), 'A', 'prg_test_units', 'pergro_longname',flag_xyfill=.true.)
        call add_default (trim(adjustl(varname)), 1, ' ')        
     enddo
  enddo

end subroutine add_fld_default_calls

!===================================================================================================
!===================================================================================================

end module physpkg
