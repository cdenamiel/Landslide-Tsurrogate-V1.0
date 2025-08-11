#!/usr/bin/python3

import numpy as np
import scipy.io
from pathlib import Path as pth

print(' ')
print("-------------------------")
print("Landslide-Tsurrogate v1.0")
print("-------------------------")
print("")
print("contact: Clea Denamiel")
print("email: clea.denamiel@live.fr")
print("")
print("Look at the User Manual and the Article Draft for more information.")
print("")
print("----------------------------------------------------------")
print("STEP 4: format the user-provided deterministic simulations")
print("----------------------------------------------------------")	
print(' ')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing this step:')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('STEP 1: Edit and run Landslide_Tsurrogate_user_input to generate the file: ../results/output_users.mat')
print('STEP 2: Run Landslide_Tsurrogate_input_parameters to generate the file: ../results/output_param.mat')
print('STEP 3: Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat')
print('        input_simus.mat: contains the matlab structure simus with zeta_max_surf[nsim,nx,ny] (maximum elevation)')
print('        ---------------  velo_max_surf[nsim,nx,ny] (maximum speed) and time_max_surf[nsim,nx,ny] (time of arrival)')
print('                         with nsim the total number of simulations corresponding to the maximum total order and')
print('                         [nx,ny] the spatial dimensions of the domain used to perform the deterministic simulations')
print('        Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat')
print('        surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models')
print('        ------------------------------  with nl the number of surrogate models to build')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('If the program and data are not edited: the program below will run for the Piton zone ')
print('                                        and the Petite Terre coastal area of the Mayotte test case.')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

def Landslide_Tsurrogate_step_4_format_simulations():	
			
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# CAUTION: the program assumes that the indices have been created in Matlab 
	#          hence index_coast = index_coast -1 and the use of flatten(order='F')
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	#-------------------------
	# Load general information
	#-------------------------
	
	# param created by Landslide_Tsurrogate_input_parameters
	my_file = pth('../results/output_param.mat')
	if my_file.is_file():
	   param_data = scipy.io.loadmat(my_file)
	   param = np.squeeze(np.array(param_data['param']))
	   maxdeg=len(param.flatten())-1
	   nsim=len(param[maxdeg]['index'][0][0].flatten())
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Run Landslide_Tsurrogate_input_param.")
	   return	
	# index_coast containing   
	my_file = pth('../data/surrogate_model_locations.mat')
	if my_file.is_file():
	   data = scipy.io.loadmat(my_file)
	   index_coast = data['index_coast'].flatten() - 1
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Produce this file as instructed above.")
	   return
	# simus_data containing the results of the deterministic simulations
	my_file = pth('../data/input_simus.mat')
	if my_file.is_file():
	   simus_data = scipy.io.loadmat(my_file)
	   zeta_max_surf = simus_data['zeta_max_surf']
	   velo_max_surf = simus_data['velo_max_surf']
	   time_max_surf = simus_data['time_max_surf']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Produce this file as instructed above.")
	   return	   
	
	#-------------------------
	# Reformat the simulations
	#-------------------------

	# Initialize the model list
	model = []
	# Loop through each simulation and extract data
	for i in range(nsim):
	    model_i = {}
	    model_i['Zeta_Z'] = zeta_max_surf[i, :, :].flatten(order='F')[index_coast]
	    model_i['Velo_Z'] = velo_max_surf[i, :, :].flatten(order='F')[index_coast]
	    model_i['Time_Zeta_Z'] = time_max_surf[i, :, :].flatten(order='F')[index_coast]
	    # Append model_i to model
	    model.append(model_i)
	
	# Save the model to a .mat file
	scipy.io.savemat('../results/output_model.mat', {'model': model})
	
	print(" ")
	print("The results have been saved in: ../results/output_model.mat")
	print(" ")
	
if __name__ == '__main__':
    Landslide_Tsurrogate_step_4_format_simulations()
