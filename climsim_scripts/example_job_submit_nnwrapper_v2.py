#!/usr/bin/env python3

# (The original python script template was provided by Walter Hannah, whannah1.)

#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
def run_cmd(cmd): print('\n'+clr.GREEN+cmd+clr.END) ; os.system(cmd); return
#---------------------------------------------------------------------------------------------------
import os, datetime, subprocess as sp, numpy as np
import shutil, glob
newcase,config,build,clean,submit,continue_run = False,False,False,False,False,False

acct = os.environ.get("MMF_NN_SLURM_ACCOUNT", "m4331")

case_prefix = 'example_job_submit_nnwrapper_v2'
# exe_refcase = ''

top_dir  = "/climsim"
case_dir = '/scratch/'
src_dir  = top_dir+'/E3SM/' 
# user_cpp = '-DMMF_ML_TRAINING' # for saving ML variables
user_cpp = '-DMMF_NN_EMULATOR' # NN hybrid test
#user_cpp = '' # do MMF

pytorch_fortran_path = '/opt/pytorch-fortran'
os.environ["pytorch_proxy_ROOT"] = pytorch_fortran_path
os.environ["pytorch_fort_proxy_ROOT"] = pytorch_fortran_path

# RESTART
runtype = 'branch' # startup, hybrid,  branch
refdate = '0002-12-30' # only works for branch (and hybrid?)
reftod = '00000' # or 21600, 43200, 64800
# clean        = True
newcase      = True
config       = True
build        = True
submit       = True
#continue_run = True
src_mod_atm  = False

debug_mode = False

dtime = 1200 # set to 0 to use a default value 

#stop_opt,stop_n,resub,walltime = 'nmonths',1, 1, '00:30:00'
stop_opt,stop_n,resub,walltime = 'nmonths',13, 0,'24:00:00'
#stop_opt,stop_n,resub,walltime = 'ndays',2, 0,'00:30:00'

ne,npg=4,2;  num_nodes=1  ; grid=f'ne{ne}pg{npg}_ne{ne}pg{npg}'
# ne,npg=30,2; num_nodes=32 ; grid=f'ne{ne}pg{npg}_ne{ne}pg{npg}'
# ne,npg=30,2; num_nodes=32 ; grid=f'ne{ne}pg{npg}_oECv3' # bi-grid for AMIP or coupled

compset,arch   = 'F2010-MMF1','GNUCPU'
# compset,arch   = 'F2010-MMF1','GNUGPU'
# compset,arch   = 'FAQP-MMF1','GNUGPU'
# compset,arch   = 'F2010-MMF1','CORI';
# (MMF1: Note that MMF_VT is tunred off for MMF_NN_EMULATOR in $E3SMROOT/components/eam/cime_config/config_component.xml)  

queue = 'regular'
#queue = 'debug'

# case_list = [case_prefix,arch,compset,grid]
case_list = [case_prefix, ]

if debug_mode: case_list.append('debug')

case='.'.join(case_list)
#---------------------------------------------------------------------------------------------------
# MMF_NN_EMULATOR
torch_model = '/storage/shared_e3sm/saved_models/wrapper/v2_mlp_wrapper.pt'
inputlength = 557
outputlength = 368
cb_nn_var_combo = 'v2'
input_rh = '.true.'
cb_spinup_step = 5
cb_strato_water_constraint = '.false.' # set .true. to use stratospheric water constraint to remove all stratospheric clouds and set dqv/dt in strato to 0
cb_partial_coupling = '.false.'
cb_do_ramp = '.false.'
cb_ramp_option = 'step'
cb_ramp_factor = 1.0
cb_ramp_step_0steps = 80
cb_ramp_step_1steps = 10

# check if MMF_ML_TRAINING is in user_cpp, then either no -DMMF_NN_EMULATOR or f_cb_partial_coupling need to be true, otherwise raise error
if 'MMF_ML_TRAINING' in user_cpp:
   if 'MMF_NN_EMULATOR' in user_cpp:
      if cb_partial_coupling == '.false.':
         raise ValueError('If MMF_NN_EMULATOR is used with MMF_ML_TRAINING, cb_partial_coupling must be true.')
#---------------------------------------------------------------------------------------------------
print('\n  case : '+case+'\n')

if 'CPU' in arch : max_mpi_per_node,atm_nthrds  =  8,1 ; max_task_per_node = 8
if 'GPU' in arch : max_mpi_per_node,atm_nthrds  =  2,8 ; max_task_per_node = 16
atm_ntasks = max_mpi_per_node*num_nodes
#---------------------------------------------------------------------------------------------------
if newcase :
   # case_scripts_dir=f'{case_dir}/{case}/case_scripts' 
   case_scripts_dir=f'{case_dir}/{case}' 
   if os.path.isdir(f'{case_dir}/{case}'): exit('\n'+clr.RED+f'This case already exists: \n{case_dir}/{case}'+clr.END+'\n')
   cmd = f'{src_dir}/cime/scripts/create_newcase -case {case} --script-root {case_scripts_dir} -compset {compset} -res {grid}  '
   if arch=='GNUCPU' : cmd += f' -mach docker-climsim -compiler gnu    -pecount {atm_ntasks}x{atm_nthrds} '
   if arch=='GNUGPU' : cmd += f' -mach docker-climsim -compiler gnugpu -pecount {atm_ntasks}x{atm_nthrds} '
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
   run_cmd('./xmlchange DOUT_S=FALSE')
   rest_option = 'ndays'
   run_cmd(f'./xmlchange REST_OPTION={rest_option}')
   run_cmd('./xmlchange REST_N=30')
   if 'max_mpi_per_node'  in locals(): run_cmd(f'./xmlchange MAX_MPITASKS_PER_NODE={max_mpi_per_node} ')
   if 'max_task_per_node' in locals(): run_cmd(f'./xmlchange MAX_TASKS_PER_NODE={max_task_per_node} ')

   # setup branch/hybrid
   if runtype == 'branch':
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_TYPE   --val {runtype}') # 'branch' won't allow change model time steps
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFDIR  --val /storage/shared_e3sm/restart_files/{refdate}-{reftod}')
      run_cmd(f'./xmlchange --file env_run.xml --id GET_REFCASE --val TRUE')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFCASE --val E3SM_ML_ne4_rerun.F2010-MMF1')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFDATE --val {refdate}')
      run_cmd(f'./xmlchange --file env_run.xml --id RUN_REFTOD  --val {reftod}')
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

&mmf_nn_emulator_nl
inputlength     = {inputlength}
outputlength    = {outputlength}
cb_nn_var_combo = '{cb_nn_var_combo}'
input_rh        = {input_rh}
cb_torch_model  = '{torch_model}'
cb_spinup_step = {cb_spinup_step}
cb_partial_coupling = {cb_partial_coupling}
cb_partial_coupling_vars = 'ptend_t', 'ptend_q0001','ptend_q0002','ptend_q0003', 'ptend_u', 'ptend_v', 'cam_out_PRECC', 'cam_out_PRECSC', 'cam_out_NETSW', 'cam_out_FLWDS', 'cam_out_SOLS', 'cam_out_SOLL', 'cam_out_SOLSD', 'cam_out_SOLLD' 
cb_do_ramp = {cb_do_ramp}
cb_ramp_option = '{cb_ramp_option}'
cb_ramp_factor = {cb_ramp_factor}
cb_ramp_step_0steps = {cb_ramp_step_0steps}
cb_ramp_step_1steps = {cb_ramp_step_1steps}
cb_strato_water_constraint = {cb_strato_water_constraint}
/

&cam_history_nl
fincl1 = 'CLDICE', 'CLDLIQ', 'DTPHYS', 'DQ1PHYS', 'DQ2PHYS', 'DQ3PHYS', 'DUPHYS'
fincl2 = 'PRECT', 'PRECC', 'FLUT', 'CLOUD', 'CLDTOT', 'CLDLOW', 'CLDMED', 'CLDHGH', 'LWCF', 'SWCF', 'LHFLX', 'SHFLX', 'TMQ', 'U850', 'T850', 'Z850', 'U500', 'T500', 'Z500', 'T', 'Q', 'U', 'V', 'CLDICE', 'CLDLIQ', 'DTPHYS', 'DQ1PHYS', 'DQ2PHYS', 'DQ3PHYS', 'DUPHYS'
avgflag_pertape = 'A','A'
nhtfrq = 0,-24
mfilt  = 0,1
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

# run_cmd(f'cp {case_dir}{exe_refcase}/build/e3sm.exe ./build/')
# run_cmd('./xmlchange BUILD_COMPLETE=TRUE')
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

