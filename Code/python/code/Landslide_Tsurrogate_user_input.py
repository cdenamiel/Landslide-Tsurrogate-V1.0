#!/usr/bin/python3

import numpy as np
import scipy.io as sio
import os

print(' ')
print("-------------------------")
print("Landslide-Tsurrogate v1.0")
print("-------------------------")
print("")
print("contact: Clea Denamiel")
print("email: clea.denamiel@live.fr")
print("")
print("For more information: ")
print(" - User Manual: ")
print(" - Publication: ")
print("")
print("-----------------------------------------------------------------------")
print("STEP 1: user-defined stochastic variables (number, type, distributions)")
print("-----------------------------------------------------------------------")
print(' ')
print(' ')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Before performing this step:')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Edit Landslide_Tsurrogate_user_input to define:')
print('                                   - the stochastic vaiables: type, number and limits of the uniform distributions')
print('                                   - the maximum total order')
print('                                   - the quadrature rule: Gauss-Patterson (GP) or Delayed Gauss-Patterson (DGP)')
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

def Landslide_Tsurrogate_user_input():	
	
	"""
	USER SPECIFICATIONS
	"""
	#---------------------------------------------------------------------------------------------
	# In this section the user defines the number of stochastic variables and their distributions.
	# In version v1.0, all the distributions are ussumed to uniform. 
	#----------------------------------------------------------------------------------------------
		
	# Number of stochastic variables used to build the surrogate models
	nmodes = 3	

	# create the vectors containing the limits of the uniform distributions
	a=np.zeros(nmodes)
	b=np.zeros(nmodes)
	# uniforme distribution of location in latitude 
	a[0] = 8584345.32
	b[0] = 8589024.15
	# uniforme distribution of volume [Mm3]
	a[1] = 1.0
	b[1] = 200.0
	# uniforme distribution of angle of friction  [degree]
	a[2] = 4.0
	b[2] = 12.0
	
	# Maximum total order of the Legendre polynomials that will be tested
	maxdeg = 6
	
	# type of nested grid option = 0 for Gauss-Patterson (GP) & option = 1 for Delayed Gauss-Patterson (DGP)
	option = 1
	
	"""
	END USER SPECIFICATIONS
	"""
	# save as mat file in order to be used later for the construction, convergence, evaluation, sensitivity of the surrogate models
	sio.savemat('../results/output_users.mat', {'nmodes': nmodes,'a': a,'b': b,'maxdeg': maxdeg,'option': option})
	
	print(" ")
	print("The results have been saved in: ../results/output_users.mat")
	print(" ")

if __name__ == '__main__':
    Landslide_Tsurrogate_user_input()
    
