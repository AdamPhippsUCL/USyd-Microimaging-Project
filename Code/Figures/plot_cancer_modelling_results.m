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
group = 'Cancer_G4';
Bools = ismember(SampleNums, eval(group));
COMP = COMP(Bools, :);

% Remove voxels with low epithelium (stroma and lumen not of interest here)
bool = (COMP(:,1)>0.4);
COMP = COMP(bool, :);


% =========== Ball+Sphere

ModelName = 'Ball+Sphere';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

% Output folder
output_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting' );

% Load parameter estimates from measured signals
measured_fs = load(fullfile(output_folder, 'Measured', ModelName, 'fs')).measured_fs;
measured_Db = load(fullfile(output_folder, 'Measured',  ModelName, 'Db')).measured_Db;
measured_R = load(fullfile(output_folder, 'Measured',  ModelName, 'R')).measured_R;

measured_fs = measured_fs(Bools);
measured_Db = measured_Db(Bools);
measured_R = measured_R(Bools);

% Load parameter estimates from predicted signals
pred_fs = load(fullfile(output_folder, 'Predicted', ModelName, 'fs')).pred_fs;
pred_Db = load(fullfile(output_folder, 'Predicted', ModelName, 'Db')).pred_Db;
pred_R = load(fullfile(output_folder, 'Predicted',  ModelName, 'R')).pred_R;

pred_fs = pred_fs(Bools);
pred_Db = pred_Db(Bools);
pred_R = pred_R(Bools);

% Only high epithelium voxels
pred_fs = pred_fs(bool);
pred_Db = pred_Db(bool);
pred_R = pred_R(bool);

measured_fs = measured_fs(bool);
measured_Db = measured_Db(bool);
measured_R = measured_R(bool);

% =========== ADC

ModelName = 'ADC';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

% Output folder
output_folder = fullfile(projectfolder, 'Outputs', 'Model Fitting' );

% Load parameter estimates from measured signals
measured_ADC = load(fullfile(output_folder, 'Measured',  ModelName, 'D')).measured_D;

measured_ADC = measured_ADC(Bools);

% Load parameter estimates from predicted signals
pred_ADC = load(fullfile(output_folder, 'Predicted', ModelName, 'D')).pred_D;

pred_ADC = pred_ADC(Bools);

% Remove voxels with low epithelium
pred_ADC = pred_ADC(bool);
measured_ADC = measured_ADC(bool);



%% SPHERE FRACTION

fs_diff = (measured_fs-pred_fs);

% Load Benign RL
RLfolder =  fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'Ball+Sphere');
BenignRL = load(fullfile(RLfolder, 'fs_BenignRL.mat')).fs_RL;
fs_bias = BenignRL(1);
fs_lowerRL = BenignRL(2);
fs_upperRL = BenignRL(3);

f=figure;
scatter(pred_fs, fs_diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on
% yline(fs_bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
yline(fs_lowerRL, '--', DisplayName='95% Limits (from benign tissue)',  color = [.1, .1, .1], LineWidth=1.2)
yline(fs_upperRL, '--', HandleVisibility="off",  color = [0.1, .1, .1], LineWidth=1.2)
legend(Location="northwest")
grid on
ylim([-0.26, 0.4])
yticks(-0.3:0.1:0.4)
xlim([-0.05, 0.35])
xticks([0:0.1:0.3])
xlabel('Predicted Sphere Fraction')
ylabel('Measured - Predicted Sphere Fraction')

ax = gca();
ax.FontSize = 12;

% Create shaded region
xlims = xlim;
xmin = xlims(1);
xmax = xlims(2);
xs = linspace(xmin, xmax, 400);

ylims = ylim;
ymin = ylims(1);
ymax = ylims(2);
ys = linspace(ymin, ymax, 400);

[X, Y] = meshgrid(xs, ys);

alphaVals = 0.1*(and(Y<fs_upperRL, Y>fs_lowerRL));

% Create base color (e.g. blue)
C = ones(size(Y,1), size(Y,2), 3);  % RGB array
C(:,:,1) = 0.1;   % red channel
C(:,:,2) = 0.1;  % green
C(:,:,3) = 0.1;  % blue (MATLAB default)

s=surf(X, Y, zeros(size(Y)), C, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'AlphaData', alphaVals, ...
    'AlphaDataMapping', 'none', ...
    HandleVisibility='off' ...
  );
uistack(s, "bottom")


f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', [group ' Residuals fs.png']))


%% BALL-COMPARTMENT DIFFUSIVITY

Db_diff = (measured_Db-pred_Db);

% Load Benign RL
RLfolder =  fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'Ball+Sphere');
BenignRL = load(fullfile(RLfolder, 'Db_BenignRL.mat')).Db_RL;
Db_bias = BenignRL(1);
Db_lowerRL = BenignRL(2);
Db_upperRL = BenignRL(3);

f=figure;
scatter(pred_Db, Db_diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on
% yline(Db_bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
yline(Db_lowerRL, '--', DisplayName='95% Limits (from benign tissue)',  color = [.1 .1 .1], LineWidth=1.2)
yline(Db_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northeast")
grid on
ylim([-0.72, 0.72])
yticks(-0.6:0.2:0.6)
xlim([0.56, 2.04])
xticks([0.6:0.2:2])
xlabel('Predicted D_{ball} (µm^2/ms)')
ylabel('Measured - Predicted D_{ball} (µm^2/ms)')

ax = gca();
ax.FontSize = 12;


% Create shaded region
xlims = xlim;
xmin = xlims(1);
xmax = xlims(2);
xs = linspace(xmin, xmax, 400);

ylims = ylim;
ymin = ylims(1);
ymax = ylims(2);
ys = linspace(ymin, ymax, 400);

[X, Y] = meshgrid(xs, ys);

alphaVals = 0.1*(and(Y<Db_upperRL, Y>Db_lowerRL));

% Create base color (e.g. blue)
C = ones(size(Y,1), size(Y,2), 3);  % RGB array
C(:,:,1) = 0.1;   % red channel
C(:,:,2) = 0.1;  % green
C(:,:,3) = 0.1;  % blue (MATLAB default)

s=surf(X, Y, zeros(size(Y)), C, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'AlphaData', alphaVals, ...
    'AlphaDataMapping', 'none', ...
    HandleVisibility='off' ...
  );
uistack(s, "bottom")


f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', [group ' Residuals Db.png']))


%% SPHERE RADIUS

R_diff = (measured_R-pred_R);

% Load Benign RL
RLfolder =  fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'Ball+Sphere');
BenignRL = load(fullfile(RLfolder, 'R_BenignRL.mat')).R_RL;
R_bias = BenignRL(1);
R_lowerRL = BenignRL(2);
R_upperRL = BenignRL(3);

f=figure;
scatter(pred_R, R_diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on

% yline(R_bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)

yline(R_lowerRL, '--', DisplayName='95% Limits (from benign tissue)',  color = [.1 .1 .1], LineWidth=1.2)
yline(R_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northeast")
grid on

ylim([-4.4, 4.6])
yticks(-6:2:6)
xlim([4, 7])
xticks([1:1:8])

xlabel('Predicted R (µm)')
ylabel('Measured - Predicted R (µm)')

ax = gca();
ax.FontSize = 12;


% Create shaded region
xlims = xlim;
xmin = xlims(1);
xmax = xlims(2);
xs = linspace(xmin, xmax, 400);

ylims = ylim;
ymin = ylims(1);
ymax = ylims(2);
ys = linspace(ymin, ymax, 400);

[X, Y] = meshgrid(xs, ys);

alphaVals = 0.1*(and(Y<R_upperRL, Y>R_lowerRL));

% Create base color (e.g. blue)
C = ones(size(Y,1), size(Y,2), 3);  % RGB array
C(:,:,1) = 0.1;   % red channel
C(:,:,2) = 0.1;  % green
C(:,:,3) = 0.1;  % blue (MATLAB default)

s=surf(X, Y, zeros(size(Y)), C, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'AlphaData', alphaVals, ...
    'AlphaDataMapping', 'none', ...
    HandleVisibility='off' ...
  );
uistack(s, "bottom")


f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', [group ' Residuals R.png']))

%% ADC

ADC_diff = (measured_ADC-pred_ADC);

% Load Benign RL
RLfolder =  fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'ADC');
BenignRL = load(fullfile(RLfolder, 'ADC_BenignRL.mat')).ADC_RL;
ADC_bias = BenignRL(1);
ADC_lowerRL = BenignRL(2);
ADC_upperRL = BenignRL(3);

f=figure;
scatter(pred_ADC, ADC_diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on
% yline(ADC_bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(ADC_lowerRL, '--', DisplayName='95% Limits (from benign tissue)',  color = [.1 .1 .1], LineWidth=1.2)
yline(ADC_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northwest")
grid on
ylim([-0.58, 0.72])
yticks(-0.8:0.2:0.8)
xlim([0.36, 1.04])
xticks([0.4:0.2:2])
xlabel('Predicted ADC (µm^2/ms)')
ylabel('Measured - Predicted ADC (µm^2/ms)')

ax = gca();
ax.FontSize = 12;

% Create shaded region
xlims = xlim;
xmin = xlims(1);
xmax = xlims(2);
xs = linspace(xmin, xmax, 400);

ylims = ylim;
ymin = ylims(1);
ymax = ylims(2);
ys = linspace(ymin, ymax, 400);

[X, Y] = meshgrid(xs, ys);

alphaVals = 0.1*(and(Y<ADC_upperRL, Y>ADC_lowerRL));%0.5*exp(-(Y-fs_bias).^2/(2*((fs_upperRL-fs_lowerRL)/3.92)^2));

% Create base color (e.g. blue)
C = ones(size(Y,1), size(Y,2), 3);  % RGB array
C(:,:,1) = 0.1;   % red channel
C(:,:,2) = 0.1;  % green
C(:,:,3) = 0.1;  % blue (MATLAB default)

s=surf(X, Y, zeros(size(Y)), C, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'AlphaData', alphaVals, ...
    'AlphaDataMapping', 'none', ...
    HandleVisibility='off' ...
  );
uistack(s, "bottom")


f.Position = [480   358   600   380];

saveas(f, fullfile(projectfolder, 'Figures', [group ' Residuals ADC.png']))