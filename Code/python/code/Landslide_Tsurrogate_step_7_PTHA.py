#!/usr/bin/python3

import scipy.io as sio
import numpy as np
from math import comb, exp
from pathlib import Path as pth
import os
import time

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
print("--------------------------------------------------------------------")
print("*******************************  PTHA  *****************************")
print("--------------------------------------------------------------------")
print(' ')
print('***************  EDIT THIS SCRIPT BEFORE RUNNING IT  ***************')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing the PTHA:')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('STEP 1: Edit and run Landslide_Tsurrogate_step_1_user_input to generate the file: ../results/output_users.mat')
print('STEP 2: Run Landslide_Tsurrogate_step_2_input_parameters to generate the file: ../results/output_param.mat')
print('STEP 3: Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat')
print('        input_simus.mat: contains zeta_max_surf[nsim,nx,ny] (maximum elevation), velo_max_surf[nsim,nx,ny] (maximum speed)')
print('        ---------------  and time_max_surf[nsim,nx,ny] (time of arrival) with nsim the total number of simulations')
print('                         corresponding to the maximum total order and [nx,ny] the spatial dimensions of the domain used to')
print('                         perform the deterministic simulations')
print('        Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat')
print('        surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models')
print('        ------------------------------  with nl the number of surrogate models to build')
print('STEP 4: Run Landslide_Tsurrogate_step_4_format_input to generate the file: ../results/output_model.mat')
print('STEP 5: Run Landslide_Tsurrogate_step_5_coefficients to generate the file: ../results/output_coeff.mat')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(' ')

def Landslide_Tsurrogate_step_7_PTHA():

        #------------
	# User Inputs
	#------------

	#--------------------------------------------------------------------------------------
	# In this section the user defines the limits of the distributions to generate the PTHA 
	#--------------------------------------------------------------------------------------

	# create the vectors containing the limits of the uniform distributions
	a0=np.zeros(3)
	b0=np.zeros(3)
	
	# chosen distances along the Piton Zone 0.0 < d < 919.0
	a0[0] = 350.0
	b0[0] = 550.0
	# chosen volume y such as 1.0 < v < 200.0
	a0[1] = 100.0
	b0[1] = 200.0
	# chosen angle of friction f such as 4.0 < f < 12.0
	a0[2] = 4.0
	b0[2] = 12.0
	
	# chosen total order maxdeg0 such as 1 < maxdeg0 < 6
	maxdeg0 = 5
	
	# chosen number of samples for the PTHA
	nw0 = 1000

        #--------------
	# Load the data
	#--------------
	
	users = pth('../results/output_users.mat')
	if users.is_file():
	   # Load the user-defined input	   
	   nmodes = np.squeeze(np.array(sio.loadmat(users)['nmodes']))	 
	   if (len(a0) != nmodes) or (len(b0) != nmodes):
	      print("The defined input for PTHA")
	      print("should be the same size than")
	      print("the number of user-defined")
	      print("stochastic variables.")
	      return
	   a = np.squeeze(np.array(sio.loadmat(users)['a']))
	   b = np.squeeze(np.array(sio.loadmat(users)['b']))
	else:
	   print("The file: ",users," does not exist.")
	   print("Edit and run Landslide_Tsurrogate_user_input.")
	   return
	my_file = pth('../results/output_coeff.mat')
	if my_file.is_file():
	   coeff = sio.loadmat(my_file)['coeff']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	   
        #--------------------------
	# Generate the PTHA
	#--------------------------
	
	start = time.time()	
	PTHA = surrogate_model_gauss_patterson_PTHA(maxdeg0, a, b, coeff, nw0, a0, b0)
	end = time.time()
	length = end - start
	
	print(" ")
	print(" ---------------------------------------------------------------")
	print("PTHA based on ",nw0," members produced in ",length, "seconds!")
	print(" ---------------------------------------------------------------")
	
	sio.savemat('../results/output_PTHA.mat',{'PTHA': PTHA})
	
	print(" ")
	print("The results have been saved in: ../results/output_PTHA.mat")
	print(" ")

def surrogate_model_gauss_patterson_PTHA(maxdeg, a, b, coeff, nw, ar, br):
        
    nmodes = len(ar)
    nx = len(coeff[0][0]['zeta'][0][0].flatten())
    
    # Build Legendre polynomials up to maxdeg (Le[0] = constant 1, Le[1] = x, etc.)
    Le = [None] * (maxdeg + 1)
    Le[0] = np.array([1.0])  # L_0 = 1
    Le[1] = np.array([1.0, 0.0])  # L_1 = x
    for n in range(2, maxdeg + 1):
        Le[n] = ((2 * n - 1) / n) * np.concatenate((Le[n - 1], [0])) - ((n - 1) / n) * np.concatenate(([0, 0], Le[n - 2]))
        
    # generate random input  
    np.random.seed(0)
    seed_01 = np.random.rand(nmodes, nw)    
    # rescaling to fit the appropriate intervals
    Zw = np.tile(ar, (nw,1)).T + np.tile((br - ar), (nw,1)).T * seed_01
               
    # Preallocate true-model matrices
    zeta_temp = np.zeros((nx, nw))
    velo_temp = np.zeros((nx, nw))
    time_temp = np.zeros((nx, nw))
    PTHA={}
    
    for alpha_norm1 in range(max(0, maxdeg - nmodes + 1), maxdeg + 1):
        # Smolyak coefficient
        C_alpha = ((-1)**(maxdeg - alpha_norm1) * comb(nmodes - 1, maxdeg - alpha_norm1))            
        # Retrieve multi-indices and PCE coefficients for this alpha_norm1
        this = coeff[0][alpha_norm1]
        alpha    = np.array(this['alpha'][0][0]) 
        zeta_hat = np.array(this['zeta'][0][0])   
        velo_hat = np.array(this['velo'][0][0])
        time_hat = np.array(this['time'][0][0])
        nww = alpha.shape[0]            
        # For each multi-index row
        for l in range(nww):
            multH = np.ones(nw)
            for n in range(nmodes):
                x_scaled = (2 * Zw[n, :] - a[n] - b[n]) / (b[n] - a[n])
                multH *= np.polyval(Le[alpha[l,n]], x_scaled)
            for i in range(nx):
                zeta_temp[i, :] += C_alpha * zeta_hat[i, l] * multH
                velo_temp[i, :] += C_alpha * velo_hat[i, l] * multH
                time_temp[i, :] += C_alpha * time_hat[i, l] * multH

    # Exponentiate to undo the log transform and get the real values of the original quantities
    PTHA['zeta'] = np.exp(zeta_temp)
    PTHA['velo'] = np.exp(velo_temp)
    PTHA['time'] = np.exp(time_temp)
    
    return PTHA

if __name__ == '__main__':
    Landslide_Tsurrogate_step_7_PTHA()
