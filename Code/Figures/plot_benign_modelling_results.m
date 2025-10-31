% Script to compare model parameters from measured and predicted signal profiles in benign tissue

clear;
projectfolder = pwd;

%% Load modelling results

% Sample groups
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };
Cancer_3 = {'4B', '4M'};
Cancer_4 = {'6N'};

folder =  fullfile(projectfolder, 'Outputs', 'Signals');
COMPOSITION = load(fullfile(folder, "COMP.mat")).COMP;
SampleNums = load(fullfile(folder, "SampleNums.mat")).SampleNums;

% Only benign samples
Bools = ismember(SampleNums, Benign);
COMP = COMPOSITION(Bools, :);
SNums = SampleNums(Bools);


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


%%  SPHERE FRACTION

fs_diff = (measured_fs-pred_fs);

% Bias
fs_bias = mean(fs_diff);

% 95% residual limits
fs_upperRL = fs_bias+1.96*std(fs_diff);
fs_lowerRL = fs_bias-1.96*std(fs_diff);

% Save residual limits
fs_RL = [fs_bias, fs_lowerRL, fs_upperRL];
RLfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'Ball+Sphere');
mkdir(RLfolder)
save(fullfile(RLfolder, 'fs_BenignRL.mat'), 'fs_RL');

f=figure;
scatter(pred_fs, fs_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP, HandleVisibility='off');
hold on
% yline(fs_bias, '-', DisplayName='Bias', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
yline(fs_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
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

f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', ['Benign Residuals Sphere Fraction.png']))


% Per sample plot

f=figure;
f.Position = [680   458   600   380];

sdiffs = [];
sindxs = [];
scolors = [];

for sindx = 1:length(Benign)

    snum = Benign{sindx};
    sbool = ismember(SNums, snum);
    sindxs = [sindxs; sindx*ones(sum(sbool), 1)];

    sdiff = fs_diff(sbool);
    scomp = COMP(sbool, :);

    sdiffs = [sdiffs; sdiff];
    scolors = [scolors; mean(scomp)];
    
end

% yline(fs_bias, '-', DisplayName='Bias', LineWidth=1.1, Alpha=0.4)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
hold on
yline(fs_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(fs_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)

% boxplot(sdiffs, sindxs,'Colors','k')

for sindx = 1:length(Benign)
    boxchart(sindxs(sindxs==sindx), sdiffs(sindxs==sindx), 'BoxFaceColor', scolors(sindx,:), 'HandleVisibility','off', MarkerStyle='.', MarkerColor=[.1, .1, .1])
    hold on
end

ylim([-0.24, 0.24])
yticks(-0.2:0.1:0.2)
xlim([0.2, 14.8])
xticks(1:20)

% Find box objects 
ax = gca();
h = findobj(ax,'Tag','Box');
set(h, 'LineWidth', 1.2)
% Loop and assign colours
for j = 1:length(h)
    set(h(j), 'Color', scolors(length(h)-j+1,:));
end

% Find outlier objects
hOut = findobj(gca,'Tag','Outliers');
set(hOut, 'MarkerEdgeColor',[0.4 0.4 0.4], 'Marker','.')

legend(Location="southeast");
ylabel('Measured - Predicted Sphere Fraction')
xlabel('Sample Number')
% title(['Sphere Fraction'])
ax.FontSize = 12;
ax.YGrid = 'on';
ax.Box = 'off';

saveas(f, fullfile(projectfolder, 'Figures', ['Samples Residuals Sphere Fraction.png']))




%% BALL-COMPARTMENT DIFFUSIVITY 

Db_diff = (measured_Db-pred_Db);

% Bias
Db_bias = mean(Db_diff);

% 95% residual limits
Db_upperRL = mean(Db_diff)+1.96*std(Db_diff);
Db_lowerRL = mean(Db_diff)-1.96*std(Db_diff);

% Save residual limits
Db_RL = [Db_bias, Db_lowerRL, Db_upperRL];
RLfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'Ball+Sphere');
mkdir(RLfolder)
save(fullfile(RLfolder, 'Db_BenignRL.mat'), 'Db_RL');

f=figure;
scatter(pred_Db, Db_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP, HandleVisibility='off');
hold on
% yline(Db_bias, '-', DisplayName='Bias', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
yline(Db_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(Db_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="southeast")
grid on
ylim([-0.69, 0.75])
yticks(-0.6:0.2:0.6)
xlim([0.56, 2.04])
xticks([0.6:0.2:2])
xlabel('Predicted D_{ball} (µm^2/ms)')
ylabel('Measured - Predicted D_{ball} (µm^2/ms)')

ax = gca();
ax.FontSize = 12;

f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', ['Benign Residuals Db.png']))


% Per sample plot

f=figure;
f.Position = [680   458   600   380];

sdiffs = [];
sindxs = [];
scolors = [];

for sindx = 1:length(Benign)

    snum = Benign{sindx};
    sbool = ismember(SNums, snum);
    sindxs = [sindxs; sindx*ones(sum(sbool), 1)];

    sdiff = Db_diff(sbool);
    scomp = COMP(sbool, :);

    sdiffs = [sdiffs; sdiff];
    scolors = [scolors; mean(scomp)];
    
end

% yline(Db_bias, '-', DisplayName='Bias', LineWidth=1.1, Alpha=0.4)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
hold on
yline(Db_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(Db_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
% boxplot(sdiffs, sindxs,'Colors','k')

for sindx = 1:length(Benign)
    boxchart(sindxs(sindxs==sindx), sdiffs(sindxs==sindx), 'BoxFaceColor', scolors(sindx,:), 'HandleVisibility','off', MarkerStyle='.', MarkerColor=[.1, .1, .1])
    hold on
end

ylim([-0.69, 0.75])
yticks(-0.6:0.2:0.6)
xlim([0.2, 14.8])
xticks(1:20)

% Find box objects 
ax = gca();
h = findobj(ax,'Tag','Box');
set(h, 'LineWidth', 1.2)
% Loop and assign colours
for j = 1:length(h)
    set(h(j), 'Color', scolors(length(h)-j+1,:));
end

% Find outlier objects
hOut = findobj(gca,'Tag','Outliers');
set(hOut, 'MarkerEdgeColor',[0.4 0.4 0.4], 'Marker','.')

legend(Location='southeast');
ylabel('Measured - Predicted D_{ball} (µm^2/ms)')
xlabel('Sample Number')
% title(['Sphere Fraction'])
ax.FontSize = 12;
ax.YGrid = 'on';
ax.Box = 'off';

saveas(f, fullfile(projectfolder, 'Figures', ['Samples Residuals Db.png']))





%% SPHERE RADIUS

% -- IGNORE VOXELS WITH HIGH FLUID CONTENT IN RESIDUAL LIMITS CALCULATION
RL_BOOL = (COMP(:,3)<0.5);


R_diff = (measured_R-pred_R);

% Bias
R_bias = mean(R_diff(RL_BOOL));

% 95% residual limits
R_upperRL = R_bias+1.96*std(R_diff(RL_BOOL));
R_lowerRL = R_bias-1.96*std(R_diff(RL_BOOL));

% Save residual limits
R_RL = [R_bias, R_lowerRL, R_upperRL];
RLfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'Ball+Sphere');
mkdir(RLfolder)
save(fullfile(RLfolder, 'R_BenignRL.mat'), 'R_RL');

f=figure;
scatter(pred_R, R_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP, HandleVisibility='off');
hold on
% yline(R_bias, '-', DisplayName='Bias', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
yline(R_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(R_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northwest")
grid on
xlim([3.8, 6.8])
xticks(2:8)
ylim([-4.9,4.2])
xlabel('Predicted R (µm)')
ylabel('Measured - Predicted R (µm)')

ax = gca();
ax.FontSize = 12;

f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', ['Benign Residuals R.png']))

% Per sample plot

f=figure;
f.Position = [680   458   600   380];

sdiffs = [];
sindxs = [];
scolors = [];

for sindx = 1:length(Benign)

    snum = Benign{sindx};
    sbool = ismember(SNums, snum);
    sindxs = [sindxs; sindx*ones(sum(sbool), 1)];

    sdiff = R_diff(sbool);
    scomp = COMP(sbool, :);

    sdiffs = [sdiffs; sdiff];
    scolors = [scolors; mean(scomp)];
    
end

% yline(R_bias, '-', DisplayName='Bias', LineWidth=1.1, Alpha=0.4)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
hold on
yline(R_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(R_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
% boxplot(sdiffs, sindxs,'Colors','k')

for sindx = 1:length(Benign)
    boxchart(sindxs(sindxs==sindx), sdiffs(sindxs==sindx), 'BoxFaceColor', scolors(sindx,:), 'HandleVisibility','off', MarkerStyle='.', MarkerColor=[.1, .1, .1])
    hold on
end

ylim([-4.9,4.2])
xlim([0.2, 14.8])
xticks(1:20)

% Find box objects 
ax = gca();
h = findobj(ax,'Tag','Box');
set(h, 'LineWidth', 1.2)
% Loop and assign colours
for j = 1:length(h)
    set(h(j), 'Color', scolors(length(h)-j+1,:));
end

% Find outlier objects
hOut = findobj(gca,'Tag','Outliers');
set(hOut, 'MarkerEdgeColor',[0.4 0.4 0.4], 'Marker','.')

lgd = legend(Location='northeast');
% lgd.Position = [0.58    0.87    0.2033    0.0118];
ylabel('Measured - Predicted R (µm)')
xlabel('Sample Number')
% title(['Sphere Fraction'])
ax.FontSize = 12;
ax.YGrid = 'on';
ax.Box = 'off';

saveas(f, fullfile(projectfolder, 'Figures', ['Samples Residuals R.png']))


%% ADC

ADC_diff = (measured_ADC-pred_ADC);

% Bias
ADC_bias = mean(ADC_diff);

% 95% residual limits
ADC_upperRL = mean(ADC_diff)+1.96*std(ADC_diff);
ADC_lowerRL = mean(ADC_diff)-1.96*std(ADC_diff);

% Save residual limits
ADC_RL = [ADC_bias, ADC_lowerRL, ADC_upperRL];
RLfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Benign RL', 'ADC');
mkdir(RLfolder)
save(fullfile(RLfolder, 'ADC_BenignRL.mat'), 'ADC_RL');

f=figure;
scatter(pred_ADC, ADC_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP, HandleVisibility='off');
hold on
% yline(ADC_bias, '-', DisplayName='Bias', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
yline(ADC_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(ADC_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="southeast")
grid on
ylim([-0.52, 0.92])
yticks(-0.8:0.2:0.8)
xlim([0.36, 2.04])
xticks([0.4:0.2:2])
xlabel('Predicted ADC (µm^2/ms)')
ylabel('Measured - Predicted ADC (µm^2/ms))')

ax = gca();
ax.FontSize = 12;

f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', ['Benign Residuals ADC.png']))


% Per sample plot

f=figure;
f.Position = [680   458   600   380];

sdiffs = [];
sindxs = [];
scolors = [];

for sindx = 1:length(Benign)

    snum = Benign{sindx};
    sbool = ismember(SNums, snum);
    sindxs = [sindxs; sindx*ones(sum(sbool), 1)];

    sdiff = ADC_diff(sbool);
    scomp = COMP(sbool, :);

    sdiffs = [sdiffs; sdiff];
    scolors = [scolors; mean(scomp)];
    
end

% yline(ADC_bias, '-', DisplayName='Bias', LineWidth=1.1, Alpha=0.4)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.1, Alpha=0.4)
hold on
yline(ADC_lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(ADC_upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
% boxplot(sdiffs, sindxs,'Colors','k')

for sindx = 1:length(Benign)
    boxchart(sindxs(sindxs==sindx), sdiffs(sindxs==sindx), 'BoxFaceColor', scolors(sindx,:), 'HandleVisibility','off', MarkerStyle='.', MarkerColor=[.1, .1, .1])
    hold on
end

xlim([0.2, 14.8])
xticks(1:20)
ylim([-0.52, 0.92])
yticks(-0.8:0.2:0.8)

% Find box objects 
ax = gca();
h = findobj(ax,'Tag','Box');
set(h, 'LineWidth', 1.2)
% Loop and assign colours
for j = 1:length(h)
    set(h(j), 'Color', scolors(length(h)-j+1,:));
end

% Find outlier objects
hOut = findobj(gca,'Tag','Outliers');
set(hOut, 'MarkerEdgeColor',[0.4 0.4 0.4], 'Marker','.')

legend(Location='southeast');
ylabel('Measured - Predicted ADC (µm^2/ms)')
xlabel('Sample Number')
% title(['Sphere Fraction'])
ax.FontSize = 12;
ax.YGrid = 'on';
ax.Box = 'off';

saveas(f, fullfile(projectfolder, 'Figures', ['Samples Residuals ADC.png']))