function Landslide_Tsurrogate_step_4_format_simulations

disp(' ')
disp('-------------------------')
disp('Landslide-Tsurrogate v1.0')
disp('-------------------------')
disp(' ')
disp('contact: Clea Denamiel')
disp('email: clea.denamiel@live.fr')
disp(' ')
disp('Look at the User Manual and the Article Draft for more information.')
disp(' ')
disp('----------------------------------------------------------')
disp('STEP 4: format the user-provided deterministic simulations')
disp('----------------------------------------------------------')
disp(' ')
disp(' ')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('Before performing this step:')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('STEP 1: Edit and run Landslide_Tsurrogate_step_1_user_input to generate the file: ../results/output_users.mat')
disp('STEP 2: Run Landslide_Tsurrogate_step_2_input_parameters to generate the file: ../results/output_param.mat')
disp('STEP 3: Run the deterministic simulations outside Landslide-Tsurrogate v1.0 and generate: ../data/input_simus.mat')
disp('        input_simus.mat: contains zeta_max_surf[nsim,nx,ny] (maximum elevation), velo_max_surf[nsim,nx,ny] (maximum speed)')
disp('        ---------------  and time_max_surf[nsim,nx,ny] (time of arrival) with nsim the total number of simulations')
disp('                         corresponding to the maximum total order and [nx,ny] the spatial dimensions of the domain used to')
disp('                         perform the deterministic simulations')
disp('        Prepare the locations where to create the surrogates and generate: ../data/surrogate_models_locations.mat')
disp('        surrogate_models_locations.mat: contains index_coast[nl] the indices where to build the surrogate models')
disp('        ------------------------------  with nl the number of surrogate models to build')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp('If the program and data are not edited: the program below will run for the Piton zone ')
disp('                                        and the Petite Terre coastal area of the Mayotte test case.')
disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
disp(' ')

%--------------------------------------------------------------------------
% Input:  
%-------
% input_simus.mat containing the results of the 
% simulations:	
%--------------------------------------------
%       zeta_max_surf:      maximum tsunami elevation 
%                           dimensions: (ns x nx x ny) => number of 
%                           simulations (ns) & spatioal dimensions (nx x ny)
%
%       velo_max_surf:      maximum tsunami speed
%                           dimensions: (ns x nx x ny)
%
%       time_max_surf:      time of arrival of the maximum tsunami elevation
%                           dimensions: (ns x nx x ny)
%
%
% surrogate_model_locations.mat containing the points where the surrogate 
% models are built: 
%--------------------------------------------
%         index_coast:      indices where to extract the surrogates
%                           dimensions: (nsm) => number of surrogate models
%             h_coast:      depths of the surrogate models
%                           dimensions: (nsm)
%           lat_coast:      latitudes of the surrogate models
%                           dimensions: (nsm)
%           lon_coast:      longitudes of the surrogate models
%                           dimensions: (nsm)
%--------------------------------------------------------------------------
% Output:
%-------
%
%               model:      structure containing the model results 
%                           dimension: (ns)
%--------------------------------------------
%              Zeta_Z:      extracted elevation corresponding to Z at stations
%                           dimensions: (nsm)
%              Velo_Z:      extracted speed corresponding to Z at stations
%                           dimensions: (nsm)
%              Time_Z:      extracted time of arrival corresponding to Z at stations
%---------------------------------------------------------------------------------------------

%--------------
% Load the data
%--------------

my_file='../results/output_param.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

my_file='../data/input_simus.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

my_file='../data/surrogate_model_locations.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

%-----------------------
% Format the simulations
%-----------------------

maxdeg = length(param); 
nsim = length(param(maxdeg).index);
for i=1:nsim
    test                 = squeeze(zeta_max_surf(i,:,:)); %#ok<*NODEF>
    model(i).Zeta_Z      = test(index_coast); %#ok<*AGROW>
    test                 = squeeze(velo_max_surf(i,:,:));
    model(i).Velo_Z      = test(index_coast);
    test                 = squeeze(time_max_surf(i,:,:));
    model(i).Time_Zeta_Z = test(index_coast);
end

%-----------------------------------------------------------------
% save as mat file in order to be used later for the construction, 
% convergence, evaluation, sensitivity of the surrogate models
%-----------------------------------------------------------------

save ../results/output_model.mat model
disp('The results have been saved in: ../results/output_model.mat')
disp(' ')

end


