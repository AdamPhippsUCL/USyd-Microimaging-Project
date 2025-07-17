% Script to plot linear model predictions and residuals

clear;
projectfolder = pwd;

%% Load scheme and measurements

SampleName = 'Multi-sample';
resultsfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement', SampleName);

% Load results
RESULTS = load(fullfile(resultsfolder, 'RESULTS.mat')).RESULTS;


%% Plot linear regression results

imgindx = 11;
bval = RESULTS(imgindx).bval;
DELTA = RESULTS(imgindx).DELTA;
R2 = RESULTS(imgindx).R2;
y = RESULTS(imgindx).y;
residuals = RESULTS(imgindx).Residuals;
X = RESULTS(1).X;

% Reformat X columns for color coding
Xnew = X;
Xnew(:,1)=X(:,2);
Xnew(:,2)=X(:,1);
X=Xnew;

f=figure;
% scatter( y-residuals, y, '.', MarkerFaceAlpha=0.25, CData=  X);
scatter( y-residuals, y,  4, 'filled', 'MarkerFaceAlpha', 0.7, CData=  X);
hold on
plot([0 0.8],[0, 0.8], color = [.1 .1 .1], LineStyle = '--', LineWidth = 1.2);
ylim([-0.02, 0.94])
xlim([-0.02, 0.94])
grid on
title(['b-value = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
ylabel('Normalized dMRI signal measurement')
xlabel('Predicted dMRI signal')
text(0.025, 0.974, ['R^2 = ' sprintf( '%0.3f', R2(1)) ' (' sprintf('%0.3f', R2(2)) ', ' sprintf('%0.3f', R2(3)) ')'], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

saveas(f, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', ['LinMod_b' num2str(bval) '_Delta' num2str(DELTA) '.png']))