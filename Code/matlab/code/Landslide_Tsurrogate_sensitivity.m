function Landslide_Tsurrogate_sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          SURROGATE MODELS FOR LANDSLIDE-GENERATED TSUNAMI
%                       CALCULATE SENSITIVITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Convergence reached at maxdeg = 5 => surrogate are built for PSA level 5
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
maxdeg = 5;
%
fprintf('*********************************************** \n')
fprintf('Sensitivity of Mayotte tsunami surrogate models \n')
fprintf('*********************************************** \n')
%
%-------------------------------
% Mayotte Nearshore Piton South 
%-------------------------------
%
fprintf('----------- \n')
fprintf('Piton South \n')
fprintf('----------- \n')
%
%----------------------------------------
% Calculate sensitivities
%----------------------------------------
load ../results/Piton_South_DGP_dim_3_order_6_coeff.mat coeff_hat
sensitivity = surrogate_model_gauss_patterson_psa_sensitivity(coeff_hat,maxdeg+1);
save ../results/Piton_South_DGP_dim_3_order_5_sensi.mat sensitivity
figure;
bar(sensitivity.ST_zeta,'stacked')
set(gcf,'Color','w')
set(gca,'FontSize',24)
title('Piton South - Elevation Sensitivity')
xlabel('Coastal Locations');
ylabel('Total sensitivity');
legend('Landslide location','Volume of the landslide','Angle of friction')
figure;
bar(sensitivity.ST_velo,'stacked')
set(gcf,'Color','w')
set(gca,'FontSize',24)
title('Piton South - Speed Sensitivity')
xlabel('Coastal Locations');
ylabel('Total sensitivity');
legend('Landslide location','Volume of the landslide','Angle of friction')
figure;
bar(sensitivity.ST_time,'stacked')
set(gcf,'Color','w')
set(gca,'FontSize',24)
title('Piton South - Time Sensitivity')
xlabel('Coastal Locations');
ylabel('Total sensitivity');
legend('Landslide location','Volume of the landslide','Angle of friction')
end
