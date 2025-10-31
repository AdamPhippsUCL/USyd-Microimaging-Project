% Display error on signal predictions (using ESL) compared to signal
% measurements (For benign samples)

clear;
projectfolder = pwd;

%% Sample and image details

% Sample groups
Cancer_3 = {'4B', '4M'};
Cancer_4 = {'6N' };
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };

group = 'Benign';

% Image
seriesindx = 11;
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
COMPOSITION = load(fullfile(folder, "COMP.mat")).COMP;
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
    case 'Cancer_3'
        Bools = ismember(SampleNums, Cancer_3);
    case 'Cancer_4'
        Bools = ismember(SampleNums, Cancer_4);
end

Pred = Predicted(Bools);
Measure = Measured(Bools);
COMP = COMPOSITION(Bools, :);


%% MAKE FIGURES

% f=figure;
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
% saveas(f, fullfile(projectfolder, 'Figures', ['Predicted vs Measured signal b' num2str(bval) '_Delta' num2str(DELTA) '.png']))
% 

%% Bland Altman

% avg = (Pred+Measure)/2;
% diff = (Measure-Pred);
% 
% % Bias
% bias = mean(diff);
% 
% % Limits of agreement
% upperLOA = mean(diff)+1.96*std(diff);
% lowerLOA = mean(diff)-1.96*std(diff);
% 
% % Save LOA
% LOA = [bias, lowerLOA, upperLOA];
% LOAfolder = fullfile(projectfolder, 'Outputs', 'Signals', SeriesDescription);
% save(fullfile(LOAfolder, 'BenignLOA.mat'), 'LOA');
% 
% f=figure;
% scatter(avg, diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP, HandleVisibility='off');
% yline(mean(diff), '-', DisplayName='Bias', LineWidth=1.2)
% hold on
% yline(mean(diff)+1.96*std(diff), '--', DisplayName='95% LOA',  color = [.1 .1 .1], LineWidth=1.2)
% yline(mean(diff)-1.96*std(diff), '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
% xlim([0, 0.8])
% ylim([-0.435, 0.435])
% legend
% grid on
% xlabel('Mean of Predicted and Measured Signal')
% ylabel('Measured Signal - Predicted Signal ')
% title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
% ax = gca();
% ax.FontSize = 12;
% f.Position = [488   242   660   400];
% 
% saveas(f, fullfile(projectfolder, 'Figures', ['Signal Bland-Altman b' num2str(bval) '_Delta' num2str(DELTA) '.png']))


%% Residuals plot

diff = (Measure-Pred);

% Bias
bias = mean(diff);

% 95% residual limits
upperRL = mean(diff)+1.96*std(diff);
lowerRL = mean(diff)-1.96*std(diff);

% Save residual limits
RL = [bias, lowerRL, upperRL];
RLfolder = fullfile(projectfolder, 'Outputs', 'Signals', SeriesDescription);
mkdir(RLfolder);
save(fullfile(RLfolder, 'BenignRL.mat'), 'RL');

f=figure;
scatter(Pred, diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=COMP, HandleVisibility='off');
hold on
% yline(bias, '-', DisplayName='Bias', LineWidth=1.2)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.2)
yline(lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
legend(Location="northwest")
grid on
ylim([-0.28, 0.28])
yticks(-0.4:0.1:0.4)
xlim([-0.02, 0.62])
xticks([0:0.1:0.6])
xlabel('Predicted Signal')
ylabel('Measured Signal - Predicted Signal')
title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
ax = gca();
ax.FontSize = 12;
f.Position = [488   242   660   400];

text(0.685, 0.95, ['R^2 = ' sprintf( '%0.3f', R2(1)) ' (' sprintf('%0.3f', R2(2)) ', ' sprintf('%0.3f', R2(3)) ')'], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

saveas(f, fullfile(projectfolder, 'Figures', ['Signal Residuals b' num2str(bval) '_Delta' num2str(DELTA) '.png']))






%% Boxplots of residual errors for each sample

f=figure;
f.Position = [488   242   660   400];

sdiffs = [];
sindxs = [];
scolors = [];

for sindx = 1:length(Benign)

    snum = Benign{sindx};
    sbool = ismember(SampleNums, snum);
    sindxs = [sindxs; sindx*ones(sum(sbool), 1)];

    spred = Predicted(sbool);
    smeasure = Measured(sbool);
    scomp = COMPOSITION(sbool, :);

    sdiffs = [sdiffs; smeasure-spred];
    scolors = [scolors; mean(scomp)];
    



end

% yline(bias, '-', DisplayName='Bias', LineWidth=1.1, Alpha=0.4)
yline(0, '-', HandleVisibility = 'off', LineWidth=1.2)

hold on
yline(lowerRL, '--', DisplayName='95% Limits',  color = [.1 .1 .1], LineWidth=1.2)
yline(upperRL, '--', HandleVisibility="off",  color = [.1 .1 .1], LineWidth=1.2)
% boxchart(sdiffs, sindxs,'Colors','k')


for sindx = 1:length(Benign)
    boxchart(sindxs(sindxs==sindx), sdiffs(sindxs==sindx), 'BoxFaceColor', scolors(sindx,:), 'HandleVisibility','off', MarkerStyle='.', MarkerColor=[.1, .1, .1])
    hold on
end

ylim([-0.28, 0.28])
yticks(-0.4:0.1:0.4)
xlim([0.2 length(Benign)+0.8])
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

legend;
ylabel('Measured Signal - Predicted Signal')
xlabel('Sample Number')
title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
ax.FontSize = 12;
ax.YGrid = 'on';
ax.Box = 'off';

saveas(f, fullfile(projectfolder, 'Figures', ['Sample Residuals b' num2str(bval) '_Delta' num2str(DELTA) '.png']))


