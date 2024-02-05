#!/usr/bin/env python3

# (The original python script template was provided by Walter Hannah, whannah1.)

#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
import os, datetime, subprocess as sp, numpy as np
import shutil, glob
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

acct = 'm4331'

case_prefix = 'FKB-rebased-test1'
# Added extra physics_state and cam_out variables.

top_dir  = os.getenv('HOME')+'/repositories'
scratch_dir = os.getenv('SCRATCH')
case_dir = scratch_dir+'/e3sm_mlt_scratch/'
src_dir  = top_dir+'/E3SM_sungdukyu/' # branch => whannah/mmf/ml-training
# user_cpp = '-DMMF_ML_TRAINING' # for saving ML variables
user_cpp = '-DCLIMSIM -DCLIMSIM_DIAG_PARTIAL -DCLIMSIMDEBUG ' # NN hybrid test
# # src_mod_atm_dir = '/global/homes/s/sungduk/repositories/ClimSim-E3SM-Hybrid/'

# RESTART
runtype = 'startup' # startup, hybrid,  branch
refdate = '0001-04-01' # only works for branch (and hybrid?)

# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
# continue_run = True
src_mod_atm  = False

debug_mode = False

dtime = 1200 # set to 0 to use a default value 

# stop_opt,stop_n,resub,walltime = 'nmonths',1, 1, '00:30:00'
stop_opt,stop_n,resub,walltime = 'ndays',15, 0,'00:30:00'

ne,npg=4,2;  num_nodes=2  ; grid=f'ne{ne}pg{npg}_ne{ne}pg{npg}'
# ne,npg=30,2; num_nodes=32 ; grid=f'ne{ne}pg{npg}_ne{ne}pg{npg}'
# ne,npg=30,2; num_nodes=32 ; grid=f'ne{ne}pg{npg}_oECv3' # bi-grid for AMIP or coupled

compset,arch   = 'F2010-MMF1','GNUGPU'
# compset,arch   = 'FAQP-MMF1','GNUGPU'
# compset,arch   = 'F2010-MMF1','CORI';
# (MMF1: Note that MMF_VT is tunred off for CLIMSIM in $E3SMROOT/components/eam/cime_config/config_component.xml)  

# queue = 'regular'
queue = 'debug'

case_list = [case_prefix,arch,compset,grid]

if debug_mode: case_list.append('debug')

case='.'.join(case_list)
#---------------------------------------------------------------------------------------------------
# CLIMSIM
f_fkb_model   = '/global/homes/s/sungduk/repositories/ClimSim-E3SM-Hybrid/model_wgts/trained_models/'\
                'backup_phase-11_retrained_models_step2_lot-152_trial_0024.best.h5.linear-out.h5.fkb.txt'
f_inp_sub     = '/global/u2/s/sungduk/repositories/ClimSim-E3SM-Hybrid/model_wgts/norm_factors/inp_sub.v2.txt'
f_inp_div     = '/global/u2/s/sungduk/repositories/ClimSim-E3SM-Hybrid/model_wgts/norm_factors/inp_div.v2.txt'
f_out_scale   = '/global/u2/s/sungduk/repositories/ClimSim-E3SM-Hybrid/model_wgts/norm_factors/out_scale.v2.txt'

#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

if 'CPU' in arch : max_mpi_per_node,atm_nthrds  = 64,1 ; max_task_per_node = 64
if 'GPU' in arch : max_mpi_per_node,atm_nthrds  =  4,8 ; max_task_per_node = 32
if arch=='CORI'  : max_mpi_per_node,atm_nthrds  = 64,1
atm_ntasks = max_mpi_per_node*num_nodes
#---------------------------------------------------------------------------------------------------
if newcase :
   # case_scripts_dir=f'{case_dir}/{case}/case_scripts' 
   case_scripts_dir=f'{case_dir}/{case}' 
   if os.path.isdir(f'{case_dir}/{case}'): exit('\n'+clr.RED+f'This case already exists: \n{case_dir}/{case}'+clr.END+'\n')
   cmd = f'{src_dir}/cime/scripts/create_newcase -case {case} --script-root {case_scripts_dir} -compset {compset} -res {grid}  '
   if arch=='GNUCPU' : cmd += f' -mach pm-cpu -compiler gnu    -pecount {atm_ntasks}x{atm_nthrds} '
   if arch=='GNUGPU' : cmd += f' -mach pm-gpu -compiler gnugpu -pecount {atm_ntasks}x{atm_nthrds} '
   if arch=='CORI'   : cmd += f' -mach cori-knl -pecount {atm_ntasks}x{atm_nthrds} '
   run_cmd(cmd)
os.chdir(f'{case_scripts_dir}')
if newcase :
   case_build_dir=f'{case_dir}/{case}/build'
   case_run_dir=f'{case_dir}/{case}/run'
   short_term_archive_root_dir=f'{case_dir}/{case}/archive'
   os.makedirs(case_build_dir, exist_ok=True)    
   os.makedirs(case_run_dir, exist_ok=True)    
   os.makedirs(short_term_archive_root_dir, exist_ok=True)    
   run_cmd(f'./xmlchange EXEROOT={case_build_dir}')
   run_cmd(f'./xmlchange RUNDIR={case_run_dir}')
   run_cmd(f'./xmlchange DOUT_S_ROOT={short_term_archive_root_dir}')
   if 'max_mpi_per_node'  in locals(): run_cmd(f'./xmlchange MAX_MPITASKS_PER_NODE={max_mpi_per_node} ')
   if 'max_task_per_node' in locals(): run_cmd(f'./xmlchange MAX_TASKS_PER_NODE={max_task_per_node} ')

   # setup branch/hybrid
   if runtype == 'branch':
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_TYPE   --val {runtype}') # 'branch' won't allow change model time steps
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFDIR  --val /pscratch/sd/s/sungduk/e3sm_mlt_scratch/REST.GNUGPU.F2010-MMF1.ne4pg2_ne4pg2/archive/rest/{refdate}-00000')
      run_cmd(f'./xmlchange --file env_run.xml --id GET_REFCASE --val TRUE')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFCASE --val REST.GNUGPU.F2010-MMF1.ne4pg2_ne4pg2')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFDATE --val {refdate}')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFTOD  --val 00000')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_STARTDATE --val {refdate}') # only used for startup or hybrid

#---------------------------------------------------------------------------------------------------
   # modify atm time step
   if dtime > 0: run_cmd(f'./xmlchange ATM_NCPL={int(24*3600/dtime)}')

   # updae user_nl_eam 
   with open("user_nl_eam", "a") as _myfile: 
       _myfile.write(f'''
&phys_ctl_nl
state_debug_checks = .true.
/

&radiation
do_aerosol_rad = .false.
/

&climsim_nl
inputlength     = 425
outputlength    = 368
cb_nn_var_combo = 'v2'
input_rh        = .false.
cb_fkb_model    = '{f_fkb_model}'
cb_inp_sub      = '{f_inp_sub}'
cb_inp_div      = '{f_inp_div}'
cb_out_scale    = '{f_out_scale}'

cb_partial_coupling = .true.
cb_partial_coupling_vars = 'ptend_t', 'ptend_q0001','ptend_q0002','ptend_q0003', 'ptend_u', 'ptend_v', 'cam_out_PRECC', 'cam_out_PRECSC', 'cam_out_NETSW', 'cam_out_FLWDS', 'cam_out_SOLS', 'cam_out_SOLL', 'cam_out_SOLSD', 'cam_out_SOLLD' 
/

&cam_history_nl
fincl2 = 'state_t_0:I:I', 'state_q0001_0:I', 'state_q0002_0:I', 'state_q0003_0:I', 'state_u_0:I', 'state_v_0:I', 'cam_out_NETSW_0:I', 'cam_out_FLWDS_0:I', 'cam_out_PRECSC_0:I', 'cam_out_PRECC_0:I', 'cam_out_SOLS_0:I', 'cam_out_SOLL_0:I', 'cam_out_SOLSD_0:I', 'cam_out_SOLLD_0:I'
fincl3 = 'state_t_1:I', 'state_q0001_1:I', 'state_q0002_1:I', 'state_q0003_1:I', 'state_u_1:I', 'state_v_1:I', 'cam_out_NETSW_1:I', 'cam_out_FLWDS_1:I', 'cam_out_PRECSC_1:I', 'cam_out_PRECC_1:I', 'cam_out_SOLS_1:I', 'cam_out_SOLL_1:I', 'cam_out_SOLSD_1:I', 'cam_out_SOLLD_1:I'
fincl4 = 'state_t_2:I', 'state_q0001_2:I', 'state_q0002_2:I', 'state_q0003_2:I', 'state_u_2:I', 'state_v_2:I', 'cam_out_NETSW_2:I', 'cam_out_FLWDS_2:I', 'cam_out_PRECSC_2:I', 'cam_out_PRECC_2:I', 'cam_out_SOLS_2:I', 'cam_out_SOLL_2:I', 'cam_out_SOLSD_2:I', 'cam_out_SOLLD_2:I'
fincl5 = 'state_t_3:I', 'state_q0001_3:I', 'state_q0002_3:I', 'state_q0003_3:I', 'state_u_3:I', 'state_v_3:I', 'cam_out_NETSW_3:I', 'cam_out_FLWDS_3:I', 'cam_out_PRECSC_3:I', 'cam_out_PRECC_3:I', 'cam_out_SOLS_3:I', 'cam_out_SOLL_3:I', 'cam_out_SOLSD_3:I', 'cam_out_SOLLD_3:I'

fincl6 = 'T:I', 'Q:I', 'CLDLIQ:I', 'CLDICE:I', 'U:I', 'V:I', 'TS:I', 'PS:I', 'LHFLX:I', 'SHFLX:I', 'SOLIN:I', 'PRECC:I', 'PRECSC:I'

nhtfrq = 0,1,1,1,1,1
mfilt  = 0,1,1,1,1,1
/

                     ''')
#---------------------------------------------------------------------------------------------------
# Copy source code modification
if src_mod_atm :
   print('Source code mods')
   dir_src_mod = f'{case_scripts_dir}/SourceMods/src.eam/'
   for kf in [*glob.glob(f'{src_mod_atm_dir}/*.F90'), *glob.glob(f'{src_mod_atm_dir}/*.xml')]:
       print(f'COPY: {kf} --> {dir_src_mod}')
       shutil.copy(kf, dir_src_mod)
#---------------------------------------------------------------------------------------------------
if config :
   cpp_defs = ''
   cpp_defs += ' '+user_cpp+' '
   if cpp_defs != '':
      run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -cppdefs \' {cpp_defs} \'  \"')
   # for ClimSim's modified namelist_definition.xml
   if src_mod_atm :
      run_cmd(f'./xmlchange --append --id CAM_CONFIG_OPTS --val \" -usr_src {dir_src_mod} \"')
   run_cmd('./xmlchange PIO_NETCDF_FORMAT=\"64bit_data\" ')
   if clean : run_cmd('./case.setup --clean')
   run_cmd('./case.setup --reset')
#---------------------------------------------------------------------------------------------------
if build : 
   if debug_mode: run_cmd('./xmlchange --file env_build.xml --id DEBUG --val TRUE ')
   if clean : run_cmd('./case.build --clean')
   run_cmd('./case.build')
#---------------------------------------------------------------------------------------------------
if submit : 
   if 'queue' in locals(): run_cmd(f'./xmlchange JOB_QUEUE={queue}')
   run_cmd(f'./xmlchange STOP_OPTION={stop_opt},STOP_N={stop_n},RESUBMIT={resub}')
   run_cmd(f'./xmlchange JOB_WALLCLOCK_TIME={walltime}')
   run_cmd(f'./xmlchange CHARGE_ACCOUNT={acct},PROJECT={acct}')
   continue_flag = 'TRUE' if continue_run else 'False'
   run_cmd(f'./xmlchange --file env_run.xml CONTINUE_RUN={continue_flag} ')   
   #-------------------------------------------------------
   run_cmd('./case.submit')
#---------------------------------------------------------------------------------------------------
# Print the case name again
print(f'\n  case : {case}\n') 
#---------------------------------------------------------------------------------------------------
#run_cmd(f'cp -v {__file__} {case_scripts_dir}/launch_script.py')
#---------------------------------------------------------------------------------------------------

