#!/usr/bin/python3

import scipy.io as sio
import numpy as np
from pathlib import Path as pth
from math import comb, exp
import os


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
print("------------------------------------------")
print("STEP 6: Evaluation of the surrogate models")
print("------------------------------------------")
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing this step:')
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

def Landslide_Tsurrogate_step_6_evaluation():

        #--------------
	# Load the data
	#--------------
	
	users = pth('../results/output_users.mat')
	if users.is_file():
	   # Load the user-defined input
	   maxdeg = sio.loadmat(users)['maxdeg'].flatten()
	   a = np.squeeze(np.array(sio.loadmat(users)['a']))
	   b = np.squeeze(np.array(sio.loadmat(users)['b']))
	else:
	   print("The file: ",users," does not exist.")
	   print("Edit and run Landslide_Tsurrogate_user_input.")
	   return
	my_file = pth('../results/output_param.mat')
	if my_file.is_file():
	   param = sio.loadmat(my_file)['param']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	my_file = pth('../results/output_model.mat')
	if my_file.is_file():	 
	   model = sio.loadmat(my_file)['model']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	my_file = pth('../results/output_coeff.mat')
	if my_file.is_file():	 
	   coeff = sio.loadmat(my_file)['coeff']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	   	   
        #---------------------------------------------
	# Generate data for evaluation and sensitivity
	#---------------------------------------------
	
	evals = surrogate_model_gauss_patterson_evals(param,model,coeff,a,b,maxdeg[0])
	sio.savemat('../results/output_evals.mat',{'evals': evals})
	
	print(" ")
	print("The results have been saved in: ../results/output_evals.mat")
	print(" ")
	
	sensi = surrogate_model_gauss_patterson_sensi(coeff,maxdeg[0]-1)
	sio.savemat('../results/output_sensi.mat',{'sensi': sensi})
	
	print(" ")
	print("The results have been saved in: ../results/output_sensi.mat")
	print(" ")

def surrogate_model_gauss_patterson_evals(param, model, coeff, a, b, maxdeg):
    """
    Evaluation of Pseudo Spectral Approximation with Gauss-Patterson Sparse grids.
    """
    
    nmodes = len(a)
    nx = len(model[0][0]['Zeta_Z'][0][0].flatten())
    
    # Build Legendre polynomials up to maxdeg (Le[0] = constant 1, Le[1] = x, etc.)
    Le = [None] * (maxdeg + 1)
    Le[0] = np.array([1.0])  # L_0 = 1
    Le[1] = np.array([1.0, 0.0])  # L_1 = x
    for n in range(2, maxdeg + 1):
        Le[n] = ((2 * n - 1) / n) * np.concatenate((Le[n - 1], [0])) - ((n - 1) / n) * np.concatenate(([0, 0], Le[n - 2]))
    # Prepare error structure as a list of dicts
    evals = []
    # Retrieve all simulations  
    Z = np.array(param[0][maxdeg]['Z'][0][0])
    index = param[0][maxdeg]['index'][0][0].flatten()
    ns = Z.shape[0]
    Z = Z.T
    # Preallocate true-model matrices
    Zeta_Z = np.empty((nx, ns))
    Velo_Z = np.empty((nx, ns))
    Time_Z = np.empty((nx, ns))
    for i in range(ns):
        entry = model[0][index[i]]
        Zeta_Z[:, i] = np.array(entry['Zeta_Z'][0][0])
        Velo_Z[:, i] = np.array(entry['Velo_Z'][0][0])
        Time_Z[:, i] = np.array(entry['Time_Zeta_Z'][0][0])
    # Store the true values as error[0]
    evals.append({
        'Zeta_Z': Zeta_Z.copy(),
        'Velo_Z': Velo_Z.copy(),
        'Time_Z': Time_Z.copy()
    })
    # Loop over increasing maximal gPCE degree
    for nmaxdeg in range(0, maxdeg+1):    
        # initialize PCE approximations
        zeta_PCE = np.zeros((nx, ns))
        velo_PCE = np.zeros((nx, ns))
        time_PCE = np.zeros((nx, ns))        
        # Smolyak summation over alpha norms
        for alpha_norm1 in range(max(0, nmaxdeg - nmodes + 1), nmaxdeg + 1):
            # Smolyak coefficient
            C_alpha = ((-1)**(nmaxdeg - alpha_norm1) 
                       * comb(nmodes-1, nmaxdeg - alpha_norm1))            
            # Retrieve multi-indices and PCE coefficients for this alpha_norm1
            this = coeff[0][alpha_norm1]
            alpha_mat = np.array(this['alpha'][0][0]) 
            zeta_hat = np.array(this['zeta'][0][0])   
            velo_hat = np.array(this['velo'][0][0])
            time_hat = np.array(this['time'][0][0])
            nw = alpha_mat.shape[0]            
            # For each multi-index row
            for l in range(nw):
                alpha_row = alpha_mat[l]                
                # Build the product of univariate polynomials
                # multH: shape (ns,)
                multH = np.ones(ns)
                for m in range(nmodes):
                    x_scaled = (2 * Z[m, :] - a[m] - b[m]) / (b[m] - a[m])
                    P = np.polyval(Le[alpha_row[m]], x_scaled)
                    multH *= P                
                # Accumulate PCE contributions
                # Note: zeta_hat[:, l] is shape (nx,), so we outer with multH
                zeta_PCE += C_alpha * np.outer(zeta_hat[:, l], multH)
                velo_PCE += C_alpha * np.outer(velo_hat[:, l], multH)
                time_PCE += C_alpha * np.outer(time_hat[:, l], multH)        
        # Exponentiate to undo the log transform and get the real values of the original quantities
        evals.append({
            'Zeta_PCE': np.exp(zeta_PCE),
            'Velo_PCE': np.exp(velo_PCE),
            'Time_PCE': np.exp(time_PCE)
        })
    
    return evals
    
def surrogate_model_gauss_patterson_sensi(coeff, maxdeg):
    
    # Dimensions
    alpha0 = np.array(coeff[0][0]['alpha'][0][0])
    nx = np.array(coeff[0][0]['zeta'][0][0]).shape[0]
    nmodes = alpha0.shape[1]
    # Preallocate total sensitivity arrays
    sensi = []
    ST_Zeta = np.full((nx, nmodes), np.nan)
    ST_Velo = np.full((nx, nmodes), np.nan)
    ST_Time = np.full((nx, nmodes), np.nan)        
    # Gather all alpha and weighted coefficients
    alpha_list = []
    coeff_zeta = []
    coeff_velo = []
    coeff_time = []
    for alpha_norm1 in range(max(0, maxdeg - nmodes + 1), maxdeg + 1):
        # Smolyak coefficient
        C_alpha = ((-1) ** (maxdeg - alpha_norm1)
                   * comb(nmodes - 1, maxdeg - alpha_norm1))        
        entry = coeff[0][alpha_norm1]
        alpha_list.append(entry['alpha'][0][0])         
        coeff_zeta.append(C_alpha * entry['zeta'][0][0])
        coeff_velo.append(C_alpha * entry['velo'][0][0])
        coeff_time.append(C_alpha * entry['time'][0][0])
    # Concatenate along the coefficient index dimension
    alpha = np.vstack(alpha_list) 
    coeff_zeta = np.hstack(coeff_zeta)
    coeff_velo = np.hstack(coeff_velo)
    coeff_time = np.hstack(coeff_time)
    # Remove purely constant terms (all zeros in alpha row)
    non_constant = ~(np.sum(alpha == 0, axis=1) == nmodes)
    alpha_nc = alpha[non_constant, :]
    zeta_nc = coeff_zeta[:, non_constant]
    velo_nc = coeff_velo[:, non_constant]
    time_nc = coeff_time[:, non_constant]
    # Total variances
    D_zeta = np.sum(zeta_nc ** 2, axis=1)
    D_velo = np.sum(velo_nc ** 2, axis=1)
    D_time = np.sum(time_nc ** 2, axis=1)
    # Compute total Sobol’ indices
    for i in range(nmodes):
        # indices of terms where alpha_nc[:, i] > 0
        indT = np.where(alpha_nc[:, i] > 0)[0]
        ST_Zeta[:, i] = np.sum(zeta_nc[:, indT] ** 2, axis=1) / D_zeta
        ST_Velo[:, i] = np.sum(velo_nc[:, indT] ** 2, axis=1) / D_velo
        ST_Time[:, i] = np.sum(time_nc[:, indT] ** 2, axis=1) / D_time
    # Store final results
    sensi.append({'ST_Zeta': ST_Zeta.copy(),
                        'ST_Velo': ST_Velo.copy(),
                        'ST_Time': ST_Time.copy()
                       })
    return sensi

if __name__ == '__main__':
    Landslide_Tsurrogate_step_6_evaluation()
    
