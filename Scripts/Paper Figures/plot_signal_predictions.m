% Script to plot linear model predictions and residuals

clear;
projectfolder = pwd;

%% Load scheme and measurements

SampleName = 'Multi-sample';
resultsfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement', SampleName);

% Load results
RESULTS = load(fullfile(resultsfolder, 'RESULTS.mat')).RESULTS;


%% Plot linear regression results

imgindx = 2;
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

figure
scatter(y, y-residuals, '*', MarkerEdgeAlpha=0.4, CData=  X);
hold on
plot(y,y, color = 'k');
ylim([-0.02, 0.82])
xlim([-0.02, 0.82])
grid on
title(['b-value = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
xlabel('Normalized dMRI signal measurement')
ylabel('Predicted dMRI signal')
text(0.03, 0.945, ['R^2 = ' sprintf( '%0.3f', R2(1)) ' (' sprintf('%0.3f', R2(2)) ', ' sprintf('%0.3f', R2(3)) ')'], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

