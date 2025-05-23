# Landslide-Tsurrogate-V1.0

**Contact**: Clea Denamiel <br>
**email (institutional)**: cdenami@irb.hr <br>
**email (permanent)**: clea.denamiel@live.fr <br>

The project contains three main folders, **Code**, **User_Manual** and **GUI**, that describe, document and share the codes of the Landslide-Tsurrogate v1.0 model.  

## Code

Contains two subfolders with the Landslide-Tsurrogate v1.0 programs and routines written in:

### matlab
### python

The workflow pipeline of the code is fully described in the associated publication (REFERENCE) and follows the steps:

**STEP 1**: Edit and run Landslide_Tsurrogate_user_input to generate the file: ../results/output_users.mat <br>
**STEP 2**: Run Landslide_Tsurrogate_input_parameters to generate the file: ../results/output_param.mat <br>
**STEP 3**: <br>
* Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat <br>
  input_simus.mat: contains zeta_max_surf[nsim,nx,ny] (maximum elevation), velo_max_surf[nsim,nx,ny] (maximum speed) and time_max_surf[nsim,nx,ny] (time of arrival) with nsim the total number of simulations corresponding to the maximum total order and [nx,ny] the spatial dimensions of the domain used to perform the deterministic simulations. <br>
*    Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat <br>
     surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models  with nl the number of surrogate models to build <br>

**STEP 4**: Run Landslide_Tsurrogate_format_input to generate the file: ../results/output_model.mat <br>
**STEP 5**: Run Landslide_Tsurrogate_psa_coefficients to generate the file: ../results/output_coeff.mat <br>
**STEP 6**: Run Landslide_Tsurrogate_surrogate_models to generate the file: ../results/output_PTHA.mat <br>

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
