#!/usr/bin/python3
import scipy.io as sio
import numpy as np
from pathlib import Path as pth
import matplotlib.pyplot as plt

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
print("-------------------------------------------")
print("*************  VISUALIZATION  *************")
print("-------------------------------------------")
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
print('STEP 6: Run Landslide_Tsurrogate_step_6_evaluation to generate the files: ../results/output_evals.mat')
print('                                                                   and ../results/output_sensi.mat')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print(' ')

def Landslide_Tsurrogate_step_6_bis_visualization():

        #--------------
	# Load the data
	#--------------
	
	users = pth('../results/output_users.mat')
	if users.is_file():
	   # Load the user-defined input
	   maxdeg = sio.loadmat(users)['maxdeg'].flatten()
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

	my_file = pth('../results/output_evals.mat')
	if my_file.is_file():	 
	   evals = sio.loadmat(my_file)['evals']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return
	
	my_file = pth('../results/output_sensi.mat')
	if my_file.is_file():	 
	   sensi = sio.loadmat(my_file)['sensi']
	else:
	   print("The file: ",my_file," does not exist.")
	   print("Generate the file as instructed above.")
	   return	   
	
	#------------------------------------
	# Convergence of the surrogate models
	#------------------------------------
		
	nsim = np.array(evals[0][0]['Zeta_Z'][0][0]).shape[0]	
	norm_err_zeta = np.full((maxdeg[0]+1, nsim), np.nan)
	norm_err_velo = np.full((maxdeg[0]+1, nsim), np.nan)
	norm_err_time = np.full((maxdeg[0]+1, nsim), np.nan)
	den_zeta = np.sum((np.array(evals[0][0]['Zeta_Z'][0][0]) - np.array(evals[0][1]['Zeta_PCE'][0][0]))**2, axis=1)
	den_velo = np.sum((np.array(evals[0][0]['Velo_Z'][0][0]) - np.array(evals[0][1]['Velo_PCE'][0][0]))**2, axis=1)
	den_time = np.sum((np.array(evals[0][0]['Time_Z'][0][0]) - np.array(evals[0][1]['Time_PCE'][0][0]))**2, axis=1)
	
	for d in range(1,maxdeg[0]+2):
	    num_zeta = np.sum((np.array(evals[0][0]['Zeta_Z'][0][0]) - np.array(evals[0][d]['Zeta_PCE'][0][0]))**2, axis=1)
	    num_velo = np.sum((np.array(evals[0][0]['Velo_Z'][0][0]) - np.array(evals[0][d]['Velo_PCE'][0][0]))**2, axis=1)
	    num_time = np.sum((np.array(evals[0][0]['Time_Z'][0][0]) - np.array(evals[0][d]['Time_PCE'][0][0]))**2, axis=1)
	    norm_err_zeta[d-1, :] = num_zeta / den_zeta
	    norm_err_velo[d-1, :] = num_velo / den_velo
	    norm_err_time[d-1, :] = num_time / den_time
	
	plot_error(norm_err_zeta.T, 'Convergence - Maximum Tsunami Elevation')
	plot_error(norm_err_velo.T, 'Convergence - Maximum Tsunami Speed')
	plot_error(norm_err_time.T, 'Convergence - Tsunami Time of Arrival')
	
	#------------------------------------
	# Accuracy of the surrogate models
	#------------------------------------	
	
	# Extract the two index arrays
	pm1=maxdeg[0]-1
	index_pm1 = param[0][maxdeg[0]-1]['index'][0][0].flatten()
	index_p = param[0][maxdeg[0]]['index'][0][0].flatten()
	# Find those in index_p that are not in index_pm1
	ind_independent = ~np.in1d(index_p, index_pm1)
	Z_true = np.array(evals[0][0]['Zeta_Z'][0][0])
	Z_pce  = np.array(evals[0][maxdeg[0]]['Zeta_PCE'][0][0])
	V_true = np.array(evals[0][0]['Velo_Z'][0][0])
	V_pce  = np.array(evals[0][maxdeg[0]]['Velo_PCE'][0][0])
	T_true  = np.array(evals[0][0]['Time_Z'][0][0])
	T_pce  = np.array(evals[0][maxdeg[0]]['Time_PCE'][0][0])
	simu_test_z = Z_true[:, ind_independent].ravel()
	pce_test_z  = Z_pce[:, ind_independent].ravel()
	simu_test_v = V_true[:, ind_independent].ravel()
	pce_test_v  = V_pce[:, ind_independent].ravel()
	simu_test_t = T_true[:, ind_independent].ravel()
	pce_test_t  = T_pce[:, ind_independent].ravel()		
	
	plot_hex(simu_test_z, pce_test_z,lim=60,vmax=15,xlabel='Determinisric simulations',ylabel='Total Order '+str(pm1),title='Accuracy - Maximum Tsunami Elevation')
	plot_hex(simu_test_v, pce_test_v,lim=50,vmax=40,xlabel='Deterministic simulations',ylabel='Total Order'+str(pm1),title='Accuracy - Maximum Tsunami Speed')
	plot_hex(simu_test_t, pce_test_t,lim=1750,vmax=200,xlabel='Deterministic simulations',ylabel='Total Order level '+str(pm1),title='Accuracy - Tsunami Time of Arrival')
	
	#------------------------------------
	# Sensitivity of the surrogate models
	#------------------------------------	
	
	ST_Zeta = np.array(sensi[0][0]['ST_Zeta'][0][0])
	ST_Velo = np.array(sensi[0][0]['ST_Velo'][0][0])
	ST_Time = np.array(sensi[0][0]['ST_Time'][0][0])
	
	locations = np.arange(ST_Zeta.shape[0])
	labels = ['Distance (Piton Zone)', 'Volume of the landslide', 'Angle of friction']
	
	plot_stacked(ST_Zeta, locations, labels, 'Sensitivity - Maximum Tsunami Elevation')
	plot_stacked(ST_Velo, locations, labels, 'Sensitivity - Maximum Tsunami Speed')
	plot_stacked(ST_Time, locations, labels, 'Sensitivity - Tsunami Time of Arrival')


def plot_error(matrix, title):
    plt.figure()
    plt.pcolormesh(matrix, shading='gouraud', cmap='Greys', vmin=0, vmax=1)
    plt.xlabel('Total Order', fontsize=14)
    plt.ylabel('Locations', fontsize=14)
    plt.title(title, fontsize=16)
    plt.colorbar(label='Normalized Squared Error')
    plt.gca().tick_params(labelsize=14)
    plt.gcf().set_facecolor('w')
    plt.show()
    
def plot_hex(simu, pce, lim, vmax, cmap='Greys', xlabel='', ylabel='', title=''):
    plt.figure()
    plt.plot([0, lim], [0, lim], 'k--', linewidth=1)
    hb = plt.hexbin(simu, pce, gridsize=100, xscale=[0, lim], yscale=[0, lim], cmap=cmap, vmin=0, vmax=vmax)
    plt.colorbar(hb, label='count')
    plt.axis('square')
    plt.xlim(0, lim)
    plt.ylim(0, lim)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.gca().tick_params(labelsize=14)
    plt.tight_layout()
    plt.show()
    
def plot_stacked(ST, locations, labels, title):
    fig, ax = plt.subplots(figsize=(10, 6))
    bottom = np.zeros(ST.shape[0])    
    # For each variable, stack its bar segment
    for i in range(ST.shape[1]):
        ax.bar(locations, ST[:, i], bottom=bottom, label=labels[i])
        bottom += ST[:, i]    
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    ax.set_title(title, fontsize=20)
    ax.set_xlabel('Coastal Locations', fontsize=16)
    ax.set_ylabel('Total sensitivity', fontsize=16)
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    plt.show()    

if __name__ == '__main__':
    Landslide_Tsurrogate_step_6_bis_visualization()
