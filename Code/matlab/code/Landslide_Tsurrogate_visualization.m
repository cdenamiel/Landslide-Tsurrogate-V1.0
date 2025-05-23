function Landslide_Tsurrogate_visualization

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

my_file='../results/output_evals.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

my_file='../results/output_sensi.mat';
if exist(my_file, 'file')
  load(my_file);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', my_file);
  uiwait(msgbox(warningMessage));
  return
end

%------------------------------------
% Convergence of the surrogate models
%------------------------------------

nsim=length(evals(1).Zeta_Z);
norm_err_zeta = nan(maxdeg+1,nsim);
norm_err_velo = nan(maxdeg+1,nsim);
norm_err_time = nan(maxdeg+1,nsim);
for dd=1:maxdeg+1
    norm_err_zeta(dd,:)      = sum((evals(1).Zeta_Z-evals(dd).Zeta_PCE).^2,2)./sum((evals(1).Zeta_Z-evals(1).Zeta_PCE).^2,2);        
    norm_err_velo(dd,:)      = sum((evals(1).Velo_Z-evals(dd).Velo_PCE).^2,2)./sum((evals(1).Velo_Z-evals(1).Velo_PCE).^2,2);
    norm_err_time(dd,:)      = sum((evals(1).Time_Z-evals(dd).Time_PCE).^2,2)./sum((evals(1).Time_Z-evals(1).Time_PCE).^2,2);
end

figure;
pcolor(norm_err_zeta');
shading interp; 
caxis([0 1]);
xlabel('PSA Level');
ylabel('Locations');
set(gca,'Fontsize',24);
set(gcf,'Color','w');
title('Convergence - Maximum Tsunami Elevation')

figure;
pcolor(norm_err_velo');
shading interp; 
caxis([0 1]);
xlabel('PSA Level');
ylabel('Locations');
set(gca,'Fontsize',24);
set(gcf,'Color','w');
title('Convergence - Maximum Tsunami Speed')

figure;
pcolor(norm_err_time');
shading interp; 
caxis([0 1]);
xlabel('PSA Level');
ylabel('Locations');
set(gca,'Fontsize',24);
set(gcf,'Color','w');
title('Convergence - Tsunami Time of Arrival')


%------------------------------------
% Accuracy of the surrogate models
%------------------------------------

index_pm1=param(maxdeg).index;
index_p=param(maxdeg+1).index;
ind_independent=(~ismember(index_p,index_pm1));

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
simu_test=evals(1).Zeta_Z(:,ind_independent);
PCE_test=evals(maxdeg).Zeta_PCE(:,ind_independent);
plot([0 60],[0 60],'k--')
hexscatter(simu_test(:),PCE_test(:),'res',100,'xlim',[0 60],'ylim',[0 60]);
axis square
ylim([0 60])
xlim([0 60])    
caxis([0 50])
xlabel('Deterministic simulations');
ylabel(['PSA level ',num2str(maxdeg-1)]);
set(gca,'Fontsize',24);
set(gcf,'Color','w');
box on
title('Accuracy - Maximum Tsunami Elevation')

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
simu_test=evals(1).Velo_Z(:,ind_independent);
PCE_test=evals(maxdeg).Velo_PCE(:,ind_independent);
plot([0 60],[0 60],'k--')
hexscatter(simu_test(:),PCE_test(:),'res',100,'xlim',[0 60],'ylim',[0 60]);
axis square
ylim([0 60])
xlim([0 60])    
caxis([0 50])
xlabel('Deterministic simulations');
ylabel(['PSA level ',num2str(maxdeg-1)]);
set(gca,'Fontsize',24);
set(gcf,'Color','w');
box on
title('Accuracy - Maximum Tsunami Speed')
%    
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
simu_test=evals(1).Time_Z(:,ind_independent);
PCE_test=evals(maxdeg).Time_PCE(:,ind_independent);
plot([0 1750],[0 1750],'k--')
hexscatter(simu_test(:),PCE_test(:),'res',100,'xlim',[0 1800],'ylim',[0 1800]);
axis square
ylim([0 1750])
xlim([0 1750])    
caxis([0 50])
xlabel('Deterministic simulations');
ylabel(['PSA level ',num2str(maxdeg-1)]);
set(gca,'Fontsize',24);
set(gcf,'Color','w');
box on
title('Accuracy - Tsunami Time of Arrival')

end
