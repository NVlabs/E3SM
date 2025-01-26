# E3SM MMF-NN-Emulator Version

## Outline

1. [Quick Start](#1-quick-start)
2. [E3SM MMF-NN-Emulator Configuration Options](#2-e3sm-mmf-nn-emulator-configuration-options)
   1. [Basic Running Modes](#21-basic-running-modes)
   2. [Partial Coupling Running Modes](#22-partical-coupling-running-modes)
3. [NN-Emulator Namelists](#3-nn-emulator-namelists)
   1. [Basic Namelists](#31-basic-namelists)
   2. [Partial Coupling Namelists](#32-nn-emulator-namelists)
4. [How to Test a Customized NN Emulator](#4-how-to-test-a-customized-nn-emulator)

---

## 1. Quick Start

To launch a hybrid E3SM MMF-NN-Emulator run, one needs to first set up the E3SM environment and download dependent files by following the instructions in the [ClimSim-Online repository](https://github.com/leap-stc/climsim-online/). These instructions will ask you to download some dependent files from a shared google drive folder [E3SM_shared](https://drive.google.com/file/d/1rH8GIx6r5rurUpzaibH3gbDksv-gViXk/view?usp=sharing) that contains some pretrained NN models and some reference E3SM-MMF simulation outputs.

Be sure to run the following command after git cloning:

```bash
git submodule update --init --recursive
```

To run a hybrid E3SM MMF-NN-Emulator simulation using our pretrained NN models used in [ClimSim-Online manuscript](https://arxiv.org/abs/2306.08754v6), you can run one of our provided example job submission scripts (```./climsim_scripts/example_job_submit_nnwrapper*.py```), e.g.,

```bash
cd climsim_scripts
python example_job_submit_v2.py
```

This example job scripts will generate 13 months of simulation data using the NN-Emulator with the pretrained MLP models that uses the "v2" input/output configuration. On NERSC you will have to run the following command before using any scripts:
```bash
module load python
```

Other two example job submission scripts (i.e., ```example_job_submit_v4.py``` and ```example_job_submit_v4_constrained.py```) use pretrained Unet models that uses the "v4" input/output configuration. Please read the section [NN-Emulator Namelists](#3-nn-emulator-namelists) for more details about 'v2' and 'v4' options.

If you run these example job script on the Perlmutter cluster, you need to modify this line `acct = os.environ.get("MMF_NN_SLURM_ACCOUNT", "m4331")` in the example job script to use your own project account instead of m4331. For personal computing resources, this line is not necessary and you don't need to modify it. The example job scripts set to create 8 tasks for the hybrid simulation. You can modify this line `if 'CPU' in arch : max_mpi_per_node,atm_nthrds  =  2,4 ; max_task_per_node = 8` in the example job script to change the number of tasks depending on your computing resources.

To use a costumized NN model, please follow section [How to Test a Customized NN Emulator](#4-how-to-test-a-customized-nn-emulator).


---

## 2. E3SM MMF-NN-Emulator configuration options

### 2.1 Basic Running Modes

Currently, there are two main CPP flags in the example job submission scripts that can control the running modes: `MMF_NN_EMULATOR`, `MMF_ML_TRAINING`. Turning them on/off will result in different model configurations as described below:

- **MMF Mode**: Turn off all the flags. Just set `user_cpp=''` in the job script (see next section).
- **MMF Mode + Saving Training Input/Output**: Turn on the `MMF_ML_TRAINING` flag. Set `user_cpp = '-DMMF_ML_TRAINING'`. This will save the necessary data every time step which we can use to construct training data for neural nets.
- **NN Mode**: Turn on `MMF_NN_EMULATOR` by set `user_cpp = '-DMMF_NN_EMULATOR'`. `MMF_NN_EMULATOR` will use NN to replace CRM calculations.

### 2.2 Partial Coupling Running Modes

- **NN Mode + Partial Coupling**: `user_cpp = '-DMMF_NN_EMULATOR'` and set namelist `cb_partial_coupling = '.true.'`. In this mode, the E3SM will do both NN and MMF calculation. You can choose to overwrite a customized set of MMF output variables by the NN outputs or a mixture of NN and MMF output (e.g. dT/dt = a * dT/dt_nn + (1-a) * dT/dt_mmf) to couple back to the rest of E3SM. In this mode, you need to set `cb_partial_coupling_vars` which will decide which variables will use NN outputs or a mixture of NN/MMF outputs. Variables not included in the `cb_partial_coupling_vars` will use the MMF output. `cb_do_ramp`, `cb_ramp_option`, and `cb_ramp_factor` can specify the customized mixture ratio of NN outputs and MMF outputs. If `cb_do_ramp=False` then we will use the `a=1` i.e. use the NN outputs to replace MMF output. If `cb_do_ramp=True` you can set `a<1` or let it depend on time. We have a few options of time schedules by setting `cb_ramp_option` and `cb_ramp_factor` see the namelist section for more details.

- **NN Mode + Partial Coupling + Saving Training Input/Output**: Similar to the previous mode but also add `-DMMF_ML_TRAINING` in the `user_cpp`. This will output input/output data which you can use for generating training data. In this mode, it will generate two sets of data. One is the GCM input/output for the (with filename containing ‘.mli.’ and ‘.mlo.’) the other is the GCM-MMF input/output (with filename containing ‘.mlis.’ and ‘.mlos.’). The difference is for the saved previous steps’ physics tendencies and advective tendencies and outputs. For the previous physics tendencies, GCM input is the prediction from the mixture of NN+MMF (may be either pure NN prediction, pure MMF prediction, or a mixture of NN + MMF prediction as specified by `cb_partial_coupling` and other related namelist as we described in the previous mode). The GCM-MMF input is the prediction from just the MMF. For the previous advective tendencies, the GCM input will be the tendencies from the `ac+dycore` while the GCM-MMF input will be the `ac+dycore+(GCM state - MMF state)/dt`. The last term is evaluated at the previous step after calling the NN/MMF module. It is nonzero only when the GCM states are not updated fully by the MMF outputs so that the GCM state is not synchronized with the CRM state. For the output, the GCM output contains the state updated by the partial coupling while the GCM-MMF output contains the state values assuming only use the MMF predictions to update. If you want to use this mode to generate training data likely the GCM-MMF output is what you want as the target because the GCM-MMF output is the pure MMF prediction. 

Note: if you use `user_cpp = '-DMMF_ML_TRAINING -DMMF_NN_EMULATOR'`, you must set `cb_partial_coupling = '.true.'` otherwise the saved `SOLIN/COSZR` will be wrong.

---

## 3. NN-Emulator namelists

### 3.1 Basic Namelists

#### torch_model
Specifies the path to the pre-trained PyTorch model.

#### inputlength
Specifies the length of NN input vector.

#### outputlength
Specifies the length of NN output vector.

#### cb_nn_var_combo
Specifies the version of input/output configuration. 'v2' is the standard input/output configuration used in the [ClimSim-Online manuscript](https://arxiv.org/abs/2306.08754v6), while 'v4' uses expanded intput features such as large-scale forcing and convection memory. Read more about these expanded input features at Section 6.3.3 in the SI of [ClimSim-Online manuscript](https://arxiv.org/abs/2306.08754v6). Currently we've only implemented 'v2' and 'v4' options. Please refer to [this file](https://github.com/leap-stc/ClimSim/blob/online_testing/climsim_utils/data_utils.py) in the [ClimSim repository](https://github.com/leap-stc/ClimSim/tree/online_testing) and look for the definition of ```self.v2_rh_inputs```, ```self.v2_outputs```, ```self.v4_inputs```, and ```self.v4_rh_outputs``` to see the detailed input/output variables for 'v2' and 'v4'. Please check the "select case" block in the ```components/eam/src/physics/cam/mmf_nn_emulator.F90``` to see how different input variables are concatenated to a single input array.

If you want to use a customized NN model, you will need to implement a new configuration in the E3SM. Please follow instructions in the section [How to Test a Customized NN Emulator](#4-how-to-test-a-customized-nn-emulator).

#### input_rh
Specifies whether to use relative humidity to replace specific humidity as input. In the 'v2' and 'v4' configurations, the 61th-120th input features are specific humidity vertical profile. If `input_rh = '.true.'`, the model will use relative humidity profile as input features. Otherwise, it will use specific humidity. To use our pretrained MLP_v2 and Unet_v4 models provided in [E3SM_shared google drive folder](https://drive.google.com/file/d/1rH8GIx6r5rurUpzaibH3gbDksv-gViXk/view?usp=sharing), you should set `input_rh = '.true.'`.

#### cb_spinup_step
It specifies how many model steps to run MMF before switching to NN. Our 'v4' configuration needs input features in previous two time steps. Therefore, in 'v4' configuration this parameter must >=3.

#### cb_strato_water_constraint
It specifies whether to use stratospheric water constraint to remove all stratospheric clouds and set dqv/dt in strato to 0. This is a microphysics constraint used in [this paper](https://arxiv.org/abs/2407.00124) that helps maintain the online stability and error growth of hybrid simulations. It is worth noting that this constraint is an approximate one. It ignores all deep overshooting stratospheric clouds which is infrequent but non-zero in actual MMF simulations.

### 3.2 Partial Coupling Namelists

#### cb_partial_coupling
It provide an option to run NN and MMF simultaneously. If it is set to '.true.', the model will run the MMF simultaneously with the NN emulator, as we described in section [Partial Coupling Running Modes](#22-partical-coupling-running-modes).

#### cb_partial_coupling_vars
If you set ```cb_partial_coupling = '.true.'```, you need to specify which variables will use NN outputs or a mixture of NN and MMF outputs. Variables not included in the ```cb_partial_coupling_vars``` will use the MMF output. For example, if you want to use NN outputs to replace the MMF output for the temperature and specific humidity tendencies, you can set ```cb_partial_coupling_vars = 'ptend_t', 'ptend_q0001'```.

#### cb_do_ramp
 If ```cb_do_ramp=False```, the model will use the NN outputs to replace MMF output in the partial coupling mode. If ```cb_do_ramp=True```, the model can use a mixture of NN and MMF outputs (e.g. dT/dt = a * dT/dt_nn + (1-a) * dT/dt_mmf) for the variables specified in ```cb_partial_coupling_vars```. you can set ```a<1``` or let it depend on time by setting ```cb_ramp_option``` and ```cb_ramp_factor```.

#### cb_ramp_option
If ```cb_do_ramp=True```, you can set ```cb_ramp_option``` to control the time schedule of the mixture ratio of NN outputs and MMF outputs. There are three options: 'constant', 'linear', 'step'. If ```cb_ramp_option = 'constant'```, ```ratio = cb_ramp_factor``` constantly in time. If ```cb_ramp_option = 'linear'```, the model will linearly ramp up the ratio value. ```ratio = cb_ramp_factor * (time/cb_ramp_linear_steps)``` when time < ```cb_ramp_linear_steps```, ```ratio = cb_ramp_factor``` when time >= ```cb_ramp_linear_steps```. Here time is counted after the initial MMF spinup steps set by ```cb_spinup_step```. If ```cb_ramp_option = 'step'```, the model will periodically switch the ratio value between ```cb_ramp_factor``` and 0.0. ```ratio = cb_ramp_factor``` for ```cb_ramp_step_1steps```, then 0.0 for ```cb_ramp_step_0steps```, then back to ```cb_ramp_factor```.

#### cb_ramp_step_1steps, cb_ramp_step_0steps
As described above, these two parameters are used to control the time schedule of the mixture ratio of NN outputs and MMF outputs when setting ```cb_ramp_option = 'step'```.

---

## 4. How to test a customized NN Emulator

In this section, we will explain how you can run a hybrid simulation with a customized NN model. 

### Using existing v2/v4 input/output configurations

You need to prepare your NN model to take input array with shape of (batch_size, inputlength) and output array with shape of (batch_size, outputlength). You need to convert your NN model into Torchscript. During training, you need to prepare your training data to be consistent with v2/v4 input/output configurations (see the definition of `self.v2_rh_inputs` and `self.v4_inputs` in [this file](https://github.com/leap-stc/ClimSim/blob/online_testing/climsim_utils/data_utils.py)).

You need to prepare your NN model to take un-normalized input features and also predict un-normalized output variables. While during model training, the input and output variables are typically normalized, you can write a wrapper model class to handle the input normalization and output un-normalization process as well as to drop unused v2 or v4 input features.

We provided instructions to reproduce our pretrained MLP_v2 and Unet_v4 models, convert them into TorchScript, and write a wrapper to deal with input/output normalization in the [ClimSim repository](https://github.com/leap-stc/ClimSim) under [this folder](https://github.com/leap-stc/ClimSim/tree/main/online_testing/model_postprocessing).

### Using customized input/output configurations

If you want to use a customized NN model with a different input configuration (you will still need to predict the same output variables with a length of 368), you need to modify `components/eam/src/physics/cam/mmf_nn_emulator.F90` to handle the new input features. In the `select case (to_lower(trim(cb_nn_var_combo)))` block, you need to add a new case option for your customized input configuration. You can refer to the existing 'v2' and 'v4' cases to see how to handle the input features.


## Author
- Zeyuan Hu, Harvard University (Previous intern at Nvidia)

## Contributors
- Akshay Subramaniam, Nvidia
- Sungduk Yu, Intel Corporation
- Walter Hannah, Lawrence Livermore National Laboratory

## References

- [ClimSim-Online: A Large Multi-scale Dataset and Framework for Hybrid ML-physics Climate Emulation](https://arxiv.org/abs/2306.08754)
- [Stable Machine-Learning Parameterization of Subgrid Processes with Real Geography and Full-physics Emulation](https://arxiv.org/abs/2407.00124)