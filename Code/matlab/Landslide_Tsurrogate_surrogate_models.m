function Landslide_Tsurrogate_surrogate_models

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
disp('--------------------------------------------------------------------')
disp('*******************************  PTHA  *****************************')
disp('--------------------------------------------------------------------')
disp(' ')
disp(' ')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Before performing this step:')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('STEP 1: Edit and run Landslide_Tsurrogate_user_input to generate the file: ../results/output_users.mat')
disp('STEP 2: Run Landslide_Tsurrogate_input_parameters to generate the file: ../results/output_param.mat')
disp('STEP 3: Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat')
disp('        input_simus.mat: contains zeta_max_surf[nsim,nx,ny] (maximum elevation), velo_max_surf[nsim,nx,ny] (maximum speed)')
disp('        ---------------  and time_max_surf[nsim,nx,ny] (time of arrival) with nsim the total number of simulations')
disp('                         corresponding to the maximum total order and [nx,ny] the spatial dimensions of the domain used to')
disp('                         perform the deterministic simulations')
disp('        Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat')
disp('        surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models')
disp('        ------------------------------  with nl the number of surrogate models to build')
disp('STEP 4: Run Landslide_Tsurrogate_format_input to generate the file: ../results/output_model.mat')
disp('STEP 5: Run Landslide_Tsurrogate_psa_coefficients to generate the file: ../results/output_coeff.mat')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INTERFACE - TO BE EDITED 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------
% In this section the user defines the limits of the distributions to generate the PTHA 
%--------------------------------------------------------------------------------------

% create the vectors containing the limits of the uniform distributions
a0 = nan(3,1);
b0 = nan(3,1);
% chosen locations in latitude y such as 8584345.32 < y < 8589024.15
a0(1) = 8589000.0;
b0(1) = 8589024.0;
% chosen volume y such as 1.0 < v < 200.0
a0(2) = 100.0;
b0(2) = 150.0;
% chosen angle of friction f such as 4.0 < f < 12.0
a0(3) = 4.0;
b0(3) = 12.0;
	
% chosen total order maxdeg0 such as 1 < maxdeg0 < 6
maxdeg0 = 5;

% chosen number of samples for the PTHA
nw0 = 20000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END USER INTERFACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------
% Load the data
%--------------

my_file='../results/output_users.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

if (length(a0) ~= nmodes) || (length(b0) ~= nmodes)
    disp('The defined input for PTHA')
    disp('should be the same size than')
	disp('the number of user-defined')
	disp('stochastic variables.')
	return
end

my_file='../results/output_coeff.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

%--------------------------
% Generate the PTHA
%--------------------------

tic
PTHA = surrogate_model_gauss_patterson_psa_PTHA(maxdeg0,a,b,coeff,nw0,a0,b0);
disp(' ')
disp(' ---------------------------------------------------------------')
disp(['PTHA based on ',num2str(nw0),' members:'])
toc
disp(' ---------------------------------------------------------------')

save ../results/output_PTHA.mat PTHA
disp(' ')
disp('The results have been saved in: ../results/output_PTHA.mat')
disp(' ')

end

function PTHA = surrogate_model_gauss_patterson_psa_PTHA(maxdeg,ar,br,coeff,nw,a,b)
% generation of the random numbers on [0 1] interval
nmodes=size(coeff(1).alpha,2);
rng(0,'twister')
seed_01 = rand(nmodes,nw);
% rescaling to fit the appropriate intervals
Zw = repmat(a,1,nw) + repmat((b-a),1,nw).*seed_01;
% Legendre Polynomials up to maxdeg
Le    = cell(maxdeg+1,1);
Le{1} = 1;       % L_1 = 1
Le{2} = [1 0];   % L_2 = x
for n = 3:maxdeg+1    
   Le{n} = ((2*n-3)/(n-1))*[Le{n-1} 0] - ((n-2)/(n-1))*[0 0 Le{n-2}];
end
% Computation of the PCE
nx=size(coeff(1).zeta,1);
zeta_temp = zeros(nx,nw);
velo_temp = zeros(nx,nw);
time_temp = zeros(nx,nw);
for alpha_norm1 = max(0,maxdeg-nmodes+1):maxdeg
    % Smolyak coefficient
    C_alpha       = (-1).^(maxdeg-alpha_norm1).*nchoosek(nmodes-1,maxdeg-alpha_norm1);
    % Retrieving the hat coefficients
    alpha         = coeff(alpha_norm1+1).alpha;
    zeta_hat      = coeff(alpha_norm1+1).zeta;
    velo_hat      = coeff(alpha_norm1+1).velo;
    time_hat      = coeff(alpha_norm1+1).time;
    nww           = size(alpha,1);
    for l = 1:nww
        multH = 1;
        for n = 1:nmodes;
            multH = multH.*polyval(Le{alpha(l,n)+1},(2*Zw(n,:)-ar(n)-br(n))./(br(n)-ar(n)));
        end
        for i = 1:nx
            zeta_temp(i,:)= squeeze(zeta_temp(i,:))+ C_alpha.*squeeze(zeta_hat(i,l))'*multH;
            velo_temp(i,:)= squeeze(velo_temp(i,:))+ C_alpha.*squeeze(velo_hat(i,l))'*multH;
            time_temp(i,:)= squeeze(time_temp(i,:))+ C_alpha.*squeeze(time_hat(i,l))'*multH;            
        end
    end
    PTHA.zeta=exp(zeta_temp);
    PTHA.velo=exp(velo_temp);
    PTHA.time=exp(time_temp);
end
end