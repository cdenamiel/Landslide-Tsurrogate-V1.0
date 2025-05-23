function Landslide_Tsurrogate_psa_evaluation

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
disp('------------------------------------------')
disp('STEP 6: Evaluation of the surrogate models')
disp('------------------------------------------')
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

my_file='../results/output_param.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

my_file='../results/output_model.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
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

%---------------------------------------------
% Generate data for evaluation and sensitivity
%---------------------------------------------

evals = surrogate_model_gauss_patterson_psa_evals(param,model,coeff,a,b,maxdeg+1);
save ../results/output_evals.mat evals
disp(' ')
disp('The results have been saved in: ../results/output_evals.mat')
disp(' ')

sensi = surrogate_model_gauss_patterson_psa_sensi(coeff,maxdeg-1);
save ../results/output_sensi.mat sensi
disp(' ')
disp('The results have been saved in: ../results/output_sensi.mat')
disp(' ')

end

function evals = surrogate_model_gauss_patterson_psa_evals(param,model,coeff,a,b,maxdeg)
%% Error of Pseudo Spectral Approximation with Gauss-Patterson Sparse grids 
% based on:
% (1) Florian Heiss, Viktor Winschel (2008): "Likelihood approximation by 
% numerical integration on sparse grids".Journal of Econometrics,
% Vol. 144, pp. 62-80.
%--------------------------------------------------------------------------
% Input:
%-------
%
% a,b:           limits of the uniform distribution such as U([a,b])
%                dimension: (nmodes)
% param:         structure containing the paramters 
%                dimension: (maxdeg+1)  
%--------------------------------------------
%         Z:     chosen samples of the stochastic variables 
%                dimensions: (nmodes) x (nwj) => dimension of the set of samples
% index:         index corresponding to the unique model run 
%
% model:         structure containing the model data 
%                dimension: (nwj)
%--------------------------------------------
%    Zeta_Z:     extracted elevation corresponding to Z at stations
%                dimensions: (nx)
% coeff:         structure containing the coefficients of the PCE 
%                dimension: (maxdeg+1)
%--------------------------------------------
%     alpha:     unique set of multi-indexes for dregree maxdeg
%                dimension: (nl) x (nmodes)
%      zeta:     polynomial coefficient for zeta at nx stations
%                dimension: (nx) x (nl) 
%--------------------------------------------------------------------------
% Output:
%--------
%
% error: structure containing the error of the PCE
%                dimension: (maxdeg+1)
%--------------------------------------------
%     Zeta_Z:    polynomial coefficient for zeta_max at nx stations
%                dimension: (nx) x (nt) x (nwj) 
%     Zeta_PCE:  polynomial coefficient for zeta_max at nx stations
%                dimension: (nx) x (nt) x (nwj) 
% -------------------------------------------------------------------------
%
nmodes      = length(a);
nx          = length(model(1).Zeta_Z);
%
% Legendre Polynomials up to maxdeg
%
Le    = cell(maxdeg,1);
Le{1} = 1;       % L_1 = 1
Le{2} = [1 0];   % L_2 = x
for n = 3:maxdeg
   Le{n} = ((2*n-3)/(n-1))*[Le{n-1} 0] - ((n-2)/(n-1))*[0 0 Le{n-2}];
end
%
% Error computation
%
evals  = struct([]);
% Retrieving all simulations
Z      = param(maxdeg).Z;
index  = param(maxdeg).index;
ns     = size(Z,1);
Z      = Z';
Zeta_Z = nan(nx,ns);
Velo_Z = nan(nx,ns);
Time_Z = nan(nx,ns);
for i = 1:ns
    Zeta_Z(:,i)  = model(index(i)).Zeta_Z;
    Velo_Z(:,i)  = model(index(i)).Velo_Z;
    Time_Z(:,i)  = model(index(i)).Time_Zeta_Z;
end
evals(1).Zeta_Z  = Zeta_Z;
evals(1).Velo_Z  = Velo_Z;
evals(1).Time_Z  = Time_Z;
for nmaxdeg = 0:maxdeg-1
    disp('*******************************')
    disp(['Calculate PCE for maxdeg = ' num2str(nmaxdeg)])
    disp('*******************************')
    % Claculate the PCE
    zeta_PCE = zeros(nx,ns);
    velo_PCE = zeros(nx,ns);
    time_PCE = zeros(nx,ns);
    for alpha_norm1 = max(0,nmaxdeg-nmodes+1):nmaxdeg
        % Smolyak coefficient
        C_alpha = (-1).^(nmaxdeg-alpha_norm1).*nchoosek(nmodes-1,nmaxdeg-alpha_norm1);
        % Retrieving the hat coefficients
        alpha    = coeff(alpha_norm1+1).alpha;
        zeta_hat = coeff(alpha_norm1+1).zeta;
        velo_hat = coeff(alpha_norm1+1).velo;
        time_hat = coeff(alpha_norm1+1).time;
        nw       = size(alpha,1);
        for l = 1:nw
            multH = 1;
            for n = 1:nmodes;
                multH = multH.*polyval(Le{alpha(l,n)+1},(2*Z(n,:)-a(n)-b(n))./(b(n)-a(n)));
            end
            for i = 1:nx        
                zeta_PCE(i,:)= squeeze(zeta_PCE(i,:))+ C_alpha.*squeeze(zeta_hat(i,l))'*multH;
                velo_PCE(i,:)= squeeze(velo_PCE(i,:))+ C_alpha.*squeeze(velo_hat(i,l))'*multH;
                time_PCE(i,:)= squeeze(time_PCE(i,:))+ C_alpha.*squeeze(time_hat(i,l))'*multH;
            end
        end
    end
    % Retrieve data for each max degree 
    evals(nmaxdeg+1).Zeta_PCE  = exp(zeta_PCE);
    evals(nmaxdeg+1).Velo_PCE  = exp(velo_PCE);
    evals(nmaxdeg+1).Time_PCE  = exp(time_PCE);
end
return
end

function sensi = surrogate_model_gauss_patterson_psa_sensi(coeff,maxdeg)
%% Sensitivity study from Pseudo Spectral Approximation method
% based on: 
% (1) Blatman & S. (2010): Sobol decomposition from PC expansions 
%--------------------------------------------------------------------------
% Input:
%-------
%
% coeff_hat: structure containing the coefficients of the PCE 
%                dimension: (maxdeg+1)
%--------------------------------------------
%     alpha:     unique set of multi-indexes for dregree maxdeg
%                dimension: (nl) x (nmodes)
%     zeta:      polynomial coefficient for zeta at nx stations
%                dimension: (nx) x (nl) 
%
%--------------------------------------------------------------------------
% Output:
%--------
%
% sensitivity: structure containing the sensitivities for Zeta
%                dimension: (maxdeg+1)
%--------------------------------------------
%     ST:        Total sensitivities at nx stations
%                dimension: (nx) x (nmodes)
% -------------------------------------------------------------------------
%
% Dimensions
nmaxdeg    = maxdeg-1;
[~,nmodes] = size(coeff(1).alpha);
[nx,~]     = size(coeff(1).zeta);
% Total Sensitivities
ST_zeta    = nan(nx,nmodes);
ST_velo    = nan(nx,nmodes);
ST_time    = nan(nx,nmodes);
%% Calculations
alpha = [];
coeff_zeta = [];
coeff_velo = [];
coeff_time = [];
for alpha_norm1 = max(0,nmaxdeg-nmodes+1):nmaxdeg
    C_alpha   = (-1).^(nmaxdeg-alpha_norm1).*nchoosek(nmodes-1,nmaxdeg-alpha_norm1);
    alpha = [alpha; coeff(alpha_norm1+1).alpha]; %#ok<AGROW>
    coeff_zeta = [coeff_zeta C_alpha.*coeff(alpha_norm1+1).zeta]; %#ok<AGROW>
    coeff_velo = [coeff_velo C_alpha.*coeff(alpha_norm1+1).velo]; %#ok<AGROW>
    coeff_time = [coeff_time C_alpha.*coeff(alpha_norm1+1).time]; %#ok<AGROW>
end
all_non_zero             = alpha(~(sum(alpha==0,2)==nmodes),:);
coef_all_non_zero_zeta   = coeff_zeta(:,~(sum(alpha==0,2)==nmodes));
D_zeta                   = sum(coef_all_non_zero_zeta.^2,2);
coef_all_non_zero_velo   = coeff_velo(:,~(sum(alpha==0,2)==nmodes));
D_velo                   = sum(coef_all_non_zero_velo.^2,2);
coef_all_non_zero_time   = coeff_time(:,~(sum(alpha==0,2)==nmodes));
D_time                   = sum(coef_all_non_zero_time.^2,2);
for i = 1:nmodes
    indT                     = find(all_non_zero(:,i) > 0);
    ST_zeta(:,i) = sum((coef_all_non_zero_zeta(:,indT)).^2,2)./D_zeta;
    ST_velo(:,i) = sum((coef_all_non_zero_velo(:,indT)).^2,2)./D_velo;    
    ST_time(:,i) = sum((coef_all_non_zero_time(:,indT)).^2,2)./D_time;
end
sensi.ST_zeta = ST_zeta;
sensi.ST_velo = ST_velo;
sensi.ST_time = ST_time;
end

