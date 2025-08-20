% Script to compare model parameters from measured and predicted signals

clear;
projectfolder = pwd;

%% Load modelling results

% Sample groups
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };
Cancer_3 = {'4B', '4M'};
Cancer_4 = {'6N'};

folder =  fullfile(projectfolder, 'Outputs', 'Signals');
COMP = load(fullfile(folder, "COMP.mat")).COMP;
SampleNums = load(fullfile(folder, "SampleNums.mat")).SampleNums;

% Only benign samples
Bools = ismember(SampleNums, Benign);
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


%% Plot results: direct comparison plot

% SPHERE FRACTION

f1=figure;
scatter(pred_fs, measured_fs,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData= COMP);
hold on
plot([0, 0.3], [0, 0.3], 'k')
grid on
xlim([-0.025,0.325])
ylim([-0.05,0.48])
xlabel('Predicted Sphere Fraction')
ylabel('Estimated Sphere Fraction')

% R2 value
SSres = sum( (pred_fs - measured_fs).^2 );
SStot = length(measured_fs)*var(measured_fs);
R2 = 1-SSres/SStot;

text(0.03, 0.945, ['R^2 = ' sprintf( '%0.3f', R2) ], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

ax = gca();
ax.FontSize = 12;

% saveas(f1, fullfile(projectfolder, 'Figures', 'Predicted vs Estimated Sphere Fraction.png'))


% BALL-COMPARTMENT DIFFUSIVITY

f2=figure;
scatter(pred_Db, measured_Db, 6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP);
hold on
plot([0.5, 2], [0.5, 2], 'k')
grid on
xlim([0.35,2.2])
ylim([0.35,2.2])
xlabel('Predicted D_{b} (x10^{-3} mm^2/s)')
ylabel('Estimated D_{b} (x10^{-3} mm^2/s)')

% R2 value
SSres = sum( (pred_Db - measured_Db).^2 );
SStot = length(measured_Db)*var(measured_Db);
R2 = 1-SSres/SStot;

text(0.03, 0.945, ['R^2 = ' sprintf( '%0.3f', R2) ], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

ax = gca();
ax.FontSize = 12;

% saveas(f2, fullfile(projectfolder, 'Figures', 'Predicted vs Estimated Db.png'))


%% Plot results: Bland-Altman

% SPHERE FRACTION

fs_avg = (pred_fs+measured_fs)/2;
fs_diff = measured_fs-pred_fs;
fs_LOA = [mean(fs_diff), mean(fs_diff)-1.96*std(fs_diff), mean(fs_diff)+1.96*std(fs_diff)];

LOA_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign LOA', ModelName);
mkdir(LOA_folder);
save( fullfile(LOA_folder,  'fs_LOA.mat'), 'fs_LOA')

f3 = figure;
scatter(fs_avg, fs_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData= COMP, HandleVisibility='off');
yline(mean(fs_diff), '-', DisplayName='Bias', LineWidth=1)
hold on
yline(mean(fs_diff)+1.96*std(fs_diff), '-.', DisplayName='95% LOA', LineWidth=1)
yline(mean(fs_diff)-1.96*std(fs_diff), '-.', HandleVisibility="off", LineWidth=1)
xlim([-0.02 0.37])
ylim([-0.35, 0.37])
xlabel('CHANGE AXES LABELS')
ylabel('CHANGE AXES LABELS')
legend
grid on
ax = gca();
ax.FontSize = 12;
% saveas(f3, fullfile(projectfolder, 'Figures', ['Bland Altman Benign Sphere Fraction.png']))


% BALL-COMPARTMENT DIFFUSIVITY

Db_avg = (pred_Db+measured_Db)/2;
Db_diff = measured_Db-pred_Db;
Db_LOA = [mean(Db_diff), mean(Db_diff)-1.96*std(Db_diff), mean(Db_diff)+1.96*std(Db_diff)];

LOA_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign LOA', ModelName);
mkdir(LOA_folder);
save( fullfile(LOA_folder,  'Db_LOA.mat'), 'Db_LOA')

f4 = figure;
scatter(Db_avg, Db_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData= COMP, HandleVisibility='off');
yline(mean(Db_diff), '-', DisplayName='Bias', LineWidth=1)
hold on
yline(mean(Db_diff)+1.96*std(Db_diff), '-.', DisplayName='95% LOA', LineWidth=1)
yline(mean(Db_diff)-1.96*std(Db_diff), '-.', HandleVisibility="off", LineWidth=1)
xlim([0.4 2.12])
ylim([-1.1, 1.1])
xlabel('CHANGE AXES LABELS')
ylabel('CHANGE AXES LABELS')
legend
grid on
ax = gca();
ax.FontSize = 12;