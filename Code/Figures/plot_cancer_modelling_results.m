% Script to compare model parameters from measured and predicted signals

clear;
projectfolder = pwd;

%% Load modelling results

% Sample groups
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };
Cancer_G3 = {'4B', '4M'};
Cancer_G4 = {'6N'};

folder =  fullfile(projectfolder, 'Outputs', 'Signals');
COMP = load(fullfile(folder, "COMP.mat")).COMP;
SampleNums = load(fullfile(folder, "SampleNums.mat")).SampleNums;

% Cancer samples
group = 'Cancer_G3';
Bools = ismember(SampleNums, eval(group));
COMP = COMP(Bools, :);

ModelName = 'Ball+Sphere';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

% Output folder
output_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting' );

% Load parameter estimates from measured signals
measured_fs = load(fullfile(output_folder, 'Measured', ModelName, 'fs')).measured_fs;
measured_Db = load(fullfile(output_folder, 'Measured',  ModelName, 'Db')).measured_Db;

measured_fs = measured_fs(Bools);
measured_Db = measured_Db(Bools);

% Load parameter estimates from predicted signals
pred_fs = load(fullfile(output_folder, 'Predicted', ModelName, 'fs')).pred_fs;
pred_Db = load(fullfile(output_folder, 'Predicted', ModelName, 'Db')).pred_Db;

pred_fs = pred_fs(Bools);
pred_Db = pred_Db(Bools);

% Remove voxels with low epithelium (stroma and lumen not of interest here)

bool = (COMP(:,1)>0.4);

pred_fs = pred_fs(bool);
measured_fs = measured_fs(bool);
pred_Db = pred_Db(bool);
measured_Db = measured_Db(bool);
COMP = COMP(bool, :);




%% Plot results: Bland-Altman
% 
% % ==== SPHERE FRACTION
% 
% fs_avg = (pred_fs+measured_fs)/2;
% fs_diff = measured_fs-pred_fs;
% 
% % Load benign LOA
% LOA_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign LOA', ModelName);
% load(fullfile(LOA_folder, 'fs_LOA.mat'));
% 
% f1 = figure;
% scatter(fs_avg, fs_diff, 14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
% yline(fs_LOA(1), '-', DisplayName='Bias (Benign)', LineWidth=1.2)
% hold on
% yline(fs_LOA(2), '--', DisplayName='95% LOA (Benign)', LineWidth=1.2)
% yline(fs_LOA(3), '--', HandleVisibility="off", LineWidth=1.2)
% xlim([-0.05 0.45])
% xticks(0:0.1:0.4)
% ylim([-0.24, 0.24])
% yticks(-0.4:0.1:0.4)
% xlabel('Mean of Measured and Predicted Sphere Fraction')
% ylabel('Measured - Predicted Sphere Fraction')
% legend(Location='northwest')
% grid on
% ax = gca();
% ax.FontSize = 12;
% f1.Position = [488   242   660   400];
% saveas(f1, fullfile(projectfolder, 'Figures', [group ' Bland-Altman Sphere Fraction.png']))
% 
% 
% 
% % ==== BALL-COMPARTMENT DIFFUSIVITY
% 
% Db_avg = (pred_Db+measured_Db)/2;
% Db_diff = measured_Db-pred_Db;
% 
% % Load Benign LOA
% LOA_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign LOA', ModelName);
% load( fullfile(LOA_folder,  'Db_LOA.mat'))
% 
% f2 = figure;
% scatter(Db_avg, Db_diff ,  14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
% yline(Db_LOA(1), '-', DisplayName='Bias (Benign)', LineWidth=1.2)
% hold on
% yline(Db_LOA(2), '--', DisplayName='95% LOA (Benign)', LineWidth=1.2)
% yline(Db_LOA(3), '--', HandleVisibility="off", LineWidth=1.2)
% xlim([0.32 2.08])
% xticks(linspace(0.2,2.2,11))
% ylim([-.72, .72])
% yticks(-0.8:0.2:0.8)
% xlabel('Mean of Measured and Predicted D_b (x10^{-3} mm^2/s)')
% ylabel('Measured - Predicted D_b (x10^{-3} mm^2/s)')
% legend(Location='northeast')
% grid on
% ax = gca();
% ax.FontSize = 12;
% 
% f2.Position = [488   242   660   400];
% saveas(f2, fullfile(projectfolder, 'Figures', [group ' Bland-Altman Db.png']))


%% Residuals plot

% SPHERE FRACTION

fs_diff = (measured_fs-pred_fs);

% Load Benign RL
RLfolder =  fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', ModelName);
BenignRL = load(fullfile(RLfolder, 'fs_BenignRL.mat')).fs_RL;
fs_bias = BenignRL(1);
fs_lowerRL = BenignRL(2);
fs_upperRL = BenignRL(3);

f=figure;
scatter(pred_fs, fs_diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on
yline(fs_bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(fs_lowerRL, '--', DisplayName='95% Residual Limits (Benign)',  color = [.1 .1 .1], LineWidth=1.2)
yline(fs_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northwest")
grid on
ylim([-0.24, 0.24])
yticks(-0.2:0.1:0.2)
xlim([-0.05, 0.35])
xticks([0:0.1:0.3])
xlabel('Predicted Sphere Fraction')
ylabel('Measured - Predicted Sphere Fraction')

ax = gca();
ax.FontSize = 12;

saveas(f, fullfile(projectfolder, 'Figures', [group ' Residuals fs.png']))


% BALL-COMPARTMENT DIFFUSIVITY

Db_diff = (measured_Db-pred_Db);

% Load Benign RL
RLfolder =  fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', ModelName);
BenignRL = load(fullfile(RLfolder, 'Db_BenignRL.mat')).Db_RL;
Db_bias = BenignRL(1);
Db_lowerRL = BenignRL(2);
Db_upperRL = BenignRL(3);

f=figure;
scatter(pred_Db, Db_diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on
yline(Db_bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(Db_lowerRL, '--', DisplayName='95% Residual Limits (Benign)',  color = [.1 .1 .1], LineWidth=1.2)
yline(Db_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northeast")
grid on
ylim([-0.72, 0.72])
yticks(-0.6:0.2:0.6)
xlim([0.56, 2.04])
xticks([0.6:0.2:2])
xlabel('Predicted Db (x10^{-3} mm^2/s)')
ylabel('Measured - Predicted Db (x10^{-3} mm^2/s)')

ax = gca();
ax.FontSize = 12;

saveas(f, fullfile(projectfolder, 'Figures', [group ' Residuals Db.png']))