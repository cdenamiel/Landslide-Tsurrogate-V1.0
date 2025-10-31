# Landslide-Tsurrogate-V1.0 -- User Manual

**Contact**: Clea Denamiel <br>
**email (institutional)**: cdenami@irb.hr <br>
**email (permanent)**: clea.denamiel@live.fr <br>

The folder contains three main sub-folders, **data**, **figures** and **results**, that describe, document and share the codes of the Landslide-Tsurrogate v1.0 model.  

## Code

Contains two subfolders with the Landslide-Tsurrogate v1.0 programs and routines written in **Matlab** & **python**.

Each of these folders conatins three subfolders: **code**, **data** and **results**. 

<img src="https://github.com/user-attachments/assets/2d911b45-22f9-42e9-9f44-0e20e79eaab4" width="300" height="300"/><br>

The code folder contains the main Landslide_Tsurrogate programs. The workflow pipeline of the code is fully described in the associated publication (article draft.pdf) and follows the steps below:

**STEP 1**: Edit and run Landslide_Tsurrogate_step_1_user_input to generate the file: ../results/output_users.mat <br>
**STEP 2**: Run Landslide_Tsurrogate_step_2_input_parameters to generate the file: ../results/output_param.mat <br>
**STEP 3**: <br>
* Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat <br>
  input_simus.mat: contains zeta_max_surf[nsim,nx,ny] (maximum elevation), velo_max_surf[nsim,nx,ny] (maximum speed) and time_max_surf[nsim,nx,ny] (time of arrival) with nsim the total number of simulations corresponding to the maximum total order and [nx,ny] the spatial dimensions of the domain used to perform the deterministic simulations. <br>
*    Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat <br>
     surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models  with nl the number of surrogate models to build <br>

**STEP 4**: Run Landslide_Tsurrogate_step_4_format_input to generate the file: ../results/output_model.mat <br>
**STEP 5**: Run Landslide_Tsurrogate_step_5_coefficients to generate the file: ../results/output_coeff.mat <br>
**STEP 6**: Run Landslide_Tsurrogate_step_6_evaluation to generate the files: ../results/output_evals.mat and ../results/output_sensi.mat  <br>
**STEP 7**: Run Landslide_Tsurrogate_step_7_PTHA to generate the file: ../results/output_PTHA.mat <br>

All these steps can be run with the provided files for the Mayotte test case: ../data/surrogate_model_locations.mat is already provided and input_simus.mat should be downloaded separately (see Code/matlab/data/README.md or Code/python/data/README.md). 

The results are created in the ../results folder.

## User Manual

The  User Manual is a Jupyter Notebook that describes, and illustrate for the Mayotte test case, the different steps and functions needed to build 
surrogate models based on gPCE for summarine landslide tsunamis.

The User Manual consists in:

- the jupyter notebook (JN): User_Manual.ipynb
- 3 different sub-folders: <br>
    * **data**: all the data used to build the surrogate models for the Mayotte test case <br>
    * **figures**: all the figures used in the JN <br>
    * **results**: an empty folder where the results from the JN are copied <br>

## GUI

Two different GUIs have been developed for the Mayotte test case using either Matlab or python.

### Matlab

Code under GUI/matlab with executables available to download for users not familiar with Matlab.

### Python  

Jupyter Widget Notebook under GUI/python with link to web application provided for users not familiar with python.
