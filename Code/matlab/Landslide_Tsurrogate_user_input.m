function Landslide_Tsurrogate_user_input

disp(' ')
disp('-------------------------')
disp('Landslide-Tsurrogate v1.0')
disp('-------------------------')
disp(' ')
disp('contact: Clea Denamiel')
disp('email: clea.denamiel@live.fr')
disp(' ')
disp('For more information: ')
disp(' - User Manual: ')
disp(' - Publication: ')
disp(' ')
disp('-----------------------------------------------------------------------')
disp('STEP 1: user-defined stochastic variables (number, type, distributions)')
disp('-----------------------------------------------------------------------')
disp(' ')
disp(' ')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Before performing this step:')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Edit Landslide_Tsurrogate_user_input to define:')
disp('                                   - the stochastic variables: type, number and limits of the uniform distributions')
disp('                                   - the maximum total order')
disp('                                   - the quadrature rule: Gauss-Patterson (GP) or Delayed Gauss-Patterson (DGP)')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INTERFACE - TO BE EDITED 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------
% Number of Stochastic Variables used to build the surrogate models
%------------------------------------------------------------------

nmodes = 3;

%------------------------------------------------------------------
% Create the vectors containing the limits of the uniform distributions
%------------------------------------------------------------------

a = nan(nmodes,1);
b = nan(nmodes,1);

% uniforme distribution of location in latitude 
a(1) = 8584345.32;
b(1) = 8589024.15;
% uniforme distribution of volume [Mm3]
a(2) = 1.0;
b(2) = 200.0;
% uniforme distribution of angle of friction  [°]
a(3) = 4.0;
b(3) = 12.0;

%-----------------------------------------------------------------
% Maximum Total Order of the Legendre polynomials that will be tested
%-----------------------------------------------------------------

maxdeg = 6;

%-----------------------------------------------------------------
% Type of nested grid option = 0 for Gauss-Patterson (GP) 
%                     option = 1 for Delayed Gauss-Patterson (DGP)
%-----------------------------------------------------------------

option = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END USER INTERFACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save ../results/output_users.mat
disp(' ')
disp('The results have been saved in: ../results/output_users.mat')
disp(' ')

end


