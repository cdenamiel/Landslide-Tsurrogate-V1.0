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

1- **Landslide_Tsurrogate_user_input** -> input saved in results/output_input.mat <br>
2- **Landslide-Tsurrogate_input_parameters** -> param saved in results/output_param.mat <br>
3- **Landslide_Tsurrogate_format_simulations** -> model saved in results/output_model.mat <br>
4- **Landslide_Tsurrogate_psa_coefficients** -> coeff saved in results/output_coeff.mat <br>
5- **Landslide_Tsurrogate_psa_evealuation** -> evals, sensi saved in results/output_evals.mat, results/output_sensi.mat <br>
   and **Landslide_Tsurrogate_visualization** -> figure produced <br>
6- **Landslide_Tsurrogate_surrogate_models** -> PTHA saved in results/output_PTHA.mat <br>

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
