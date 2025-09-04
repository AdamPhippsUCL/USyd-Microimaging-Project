% Display error on signal predictions (using ESL) compared to signal
% measurements (For benign samples)

clear;
projectfolder = pwd;

%% Sample and image details

% Sample groups
Cancer_3 = {'4B', '4M'};
Cancer_4 = {'6N' };
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };

group = 'Cancer_G4';

% Image
seriesindx =11;
SeriesDescriptions = {
    'SE_b0_SPOIL5% (DS)',...
    'STEAM_ShortDELTA_15 (DS)',...
    'STEAM_ShortDELTA_20 (DS)',...
    'STEAM_ShortDELTA_30 (DS)',...
    'STEAM_ShortDELTA_40 (DS)',...
    'STEAM_ShortDELTA_50 (DS)',...
    'STEAM_LongDELTA_40 (DS)',...
    'STEAM_LongDELTA_60 (DS)',...
    'STEAM_LongDELTA_80 (DS)',...
    'STEAM_LongDELTA_100 (DS)',...
    'STEAM_LongDELTA_120 (DS)'...
};
SeriesDescription = SeriesDescriptions{seriesindx};

scheme = load(fullfile(projectfolder, "Schemes", "20250224_UQ4 AllDELTA.mat")).scheme;
bval = scheme(seriesindx).bval;
DELTA = scheme(seriesindx).DELTA;

% Load signals + extras
folder =  fullfile(projectfolder, 'Outputs', 'Signals');
COMP = load(fullfile(folder, "COMP.mat")).COMP;
SampleNums = load(fullfile(folder, "SampleNums.mat")).SampleNums;
Measured = load(fullfile(folder, SeriesDescription, "Measured.mat")).Measured;
Predicted = load(fullfile(folder, SeriesDescription, "Predicted.mat")).Predicted;

% R2 value
RESULTS = load(fullfile(projectfolder, 'Outputs', 'ESL signal estimation', 'Multi-sample', 'RESULTS.mat')).RESULTS;
R2 = RESULTS(seriesindx).R2;
clear RESULTS

% For samples in group only
switch group
    case 'Benign'
        Bools = ismember(SampleNums, Benign);
    case 'Cancer_G3'
        Bools = ismember(SampleNums, Cancer_3);
    case 'Cancer_G4'
        Bools = ismember(SampleNums, Cancer_4);
end

Pred = Predicted(Bools);
Measure = Measured(Bools);
COMP = COMP(Bools, :);

% Remove voxels with low epithelium (stroma and lumen not of interest here)
bool = (COMP(:,1)>0.4);
Pred = Pred(bool);
Measure = Measure(bool);
COMP = COMP(bool, :);

% %% MAKE FIGURES
% 
% f=figure;
% 
% scatter(Pred, Measure,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData= COMP);
% hold on
% plot([0 0.8],[0, 0.8], color = [.1 .1 .1], LineStyle = '--', LineWidth = 1.2);
% 
% ylim([-0.02, 0.8])
% xlim([-0.02, 0.8])
% grid on
% xlabel('Predicted Signal')
% ylabel('Measured Signal')
% title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
% 
% text(0.025, 0.97, ['R^2 = ' sprintf( '%0.3f', R2(1)) ' (' sprintf('%0.3f', R2(2)) ', ' sprintf('%0.3f', R2(3)) ')'], ...
%     'Units', 'normalized', ...
%     'VerticalAlignment', 'top', ...
%     'HorizontalAlignment', 'left', ...
%     'BackgroundColor', 'white', ...
%     'EdgeColor', 'black');  % Optional border
% 
% f.Position = [488   242   660   400];
% ax = gca();
% ax.FontSize = 12;
% 
% % saveas(f, fullfile(projectfolder, 'Figures', ['Predicted vs Measured signal b' num2str(bval) '_Delta' num2str(DELTA) '.png']))


% %% Bland Altman
% 
% avg = (Pred+Measure)/2;
% diff = (Measure-Pred);
% 
% % Load Benign LOA
% LOAfolder = fullfile(projectfolder, 'Outputs', 'Signals', SeriesDescription);
% BenignLOA = load(fullfile(LOAfolder, 'BenignLOA.mat')).LOA;
% 
% f=figure;
% scatter(avg, diff ,  14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
% yline(BenignLOA(1), '-', DisplayName='Bias (Benign)', LineWidth=1.2)
% hold on
% yline(BenignLOA(2), '--', DisplayName='95% LOA (Benign)', LineWidth=1.2)
% yline(BenignLOA(3), '--', HandleVisibility="off", LineWidth=1.2)
% xlim([0, 0.8])
% ylim([-0.435, 0.435])
% legend(Location="northwest")
% grid on
% xlabel('Mean of Predicted and Measured Signal')
% ylabel('Measured Signal - Predicted Signal ')
% title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
% ax = gca();
% ax.FontSize = 12;
% f.Position = [488   242   660   400];
% saveas(f, fullfile(projectfolder, 'Figures', [group ' Signal Bland-Altman b' num2str(bval) '_Delta' num2str(DELTA) '.png']))
% 


%% Residuals plot

diff = (Measure-Pred);

% Load Benign RL
RLfolder = fullfile(projectfolder, 'Outputs', 'Signals', SeriesDescription);
BenignRL = load(fullfile(RLfolder, 'BenignRL.mat')).RL;
bias = BenignRL(1);
lowerRL = BenignRL(2);
upperRL = BenignRL(3);


f=figure;
scatter(Pred, diff ,   14, 'filled', 'MarkerFaceAlpha', 1, CData=COMP, HandleVisibility='off')
hold on
yline(bias, '-', DisplayName='Bias (Benign)', LineWidth=1.2)
yline(lowerRL, '--', DisplayName='95% Residual Limits (Benign)',  color = [.1 .1 .1], LineWidth=1.2)
yline(upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northwest")
grid on
ylim([-0.42, 0.42])
yticks(-0.4:0.1:0.4)
xlim([-0.05, 0.55])
xticks([0:0.1:0.5])
xlabel('Predicted Signal')
ylabel('Measured Signal - Predicted Signal')
title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
ax = gca();
ax.FontSize = 12;


f.Position = [680   458   600   380];

saveas(f, fullfile(projectfolder, 'Figures', [group ' Signal Residuals b' num2str(bval) '_Delta' num2str(DELTA) '.png']))

