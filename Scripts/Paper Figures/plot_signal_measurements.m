% Script to plot signal measurements

clear;
projectfolder = pwd;


%% Load scheme and measurements

SampleName = 'Multi-sample';
resultsfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement', SampleName);

signals = load(fullfile(resultsfolder, 'signals.mat')).signals;
scheme = load(fullfile(resultsfolder, 'scheme.mat')).scheme;


%% Create figure

f1=figure;
bshift=5;

% Short Delta

indices = 2:6;
s = signals(:,indices,:);
bvals = [scheme(indices).bval];
lw = 0.8;
% Display error bars with 95% confidence intervals from botstrapping
errorbar(bvals-1*bshift, s(2,:,1), s(2,:,1)-s(2,:,3), s(2,:,4)-s(2,:,1), '-*', LineWidth = lw, color=[0.8500 0.3250 0.0980], DisplayName='E (Short \Delta)');
hold on
errorbar(bvals-1*bshift, s(1,:,1), s(1,:,1)-s(1,:,3), s(1,:,4)-s(1,:,1), '-*', LineWidth = lw, color=[0.4660 0.6740 0.1880], DisplayName='S (Short \Delta)');
errorbar(bvals-1*bshift, s(3,:,1), s(3,:,1)-s(3,:,3), s(3,:,4)-s(3,:,1),'-*', LineWidth = lw, color=[0 0.4470 0.7410], DisplayName = 'L (Short \Delta)');

% Long Delta
indices = 7:11;
s = signals(:,indices,:);

errorbar(bvals+1*bshift, s(2,:,1), s(2,:,1)-s(2,:,3), s(2,:,4)-s(2,:,1),  '-.*', LineWidth = lw, color=[0.8500 0.3250 0.0980], DisplayName='E (Long \Delta)');
errorbar(bvals+1*bshift, s(1,:,1), s(1,:,1)-s(1,:,3), s(1,:,4)-s(1,:,1), '-.*', LineWidth = lw, color=[0.4660 0.6740 0.1880], DisplayName='S (Long \Delta)');
errorbar(bvals+1*bshift, s(3,:,1), s(3,:,1)-s(3,:,3), s(3,:,4)-s(3,:,1), '-.*', LineWidth = lw, color=[0 0.4470 0.7410], DisplayName = 'L (Long \Delta)');

xticks(bvals); 
xticklabels(bvals)
% ylim([0.18,0.68])
ylim([-0.05, 0.75])
ylabel('Normalized dMRI signal measurement')
yticks(linspace(0.0, 1, 11))
xlim([800,2200])
xlabel('b-value (s/mm^{2})')
xticks(bvals)
xticklabels(["1000", "1250", "1500", "1750", "2000"])
grid on
legend('NumColumns', 2);
ax=gca();
ax.FontSize=12;
f1.Position = [488   242   720   600];
saveas(f1, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', 'Signal_Measurements.png'))

%% ES contrast

f2=figure;
f2.Position = [680   458   720   420];
bshift = 15;
% Short Delta
indices = 2:6;
s = signals(:,indices,:);
bvals = [scheme(indices).bval];
% errorbar(bvals-1*bshift, s(2,:,1)-s(1,:,1), s(2,:,1)-s(2,:,3) + s(1,:,1)-s(1,:,3), s(2,:,4)-s(2,:,1) + s(1,:,4)-s(1,:,1), '--*', DisplayName = 'Short \Delta')
short_contrast = s(2,:,1)-s(1,:,1);
short_SE = sqrt(s(2,:,2).^2 + s(1,:,2).^2);
short_CI = 1.96*short_SE;

errorbar(bvals-1*bshift, short_contrast, short_CI, '-*', LineWidth = lw, color = 'k' , DisplayName = 'Short \Delta')
hold on


% Long Delta
indices = 7:11;
s = signals(:,indices,:);
% errorbar(bvals+1*bshift, s(2,:,1)-s(1,:,1), s(2,:,1)-s(2,:,3) + s(1,:,1)-s(1,:,3), s(2,:,4)-s(2,:,1) + s(1,:,4)-s(1,:,1), '--*', DisplayName = 'Long \Delta')

long_contrast = s(2,:,1)-s(1,:,1);
long_SE = sqrt(s(2,:,2).^2 + s(1,:,2).^2);
long_CI = 1.96*long_SE;

errorbar(bvals+1*bshift, long_contrast, long_CI, '--*', LineWidth = lw, color = 'k' ,DisplayName = 'Long \Delta')
% 95% confidence intervals calculated unsing standard error of signal
% measurements


ylim([0.02, 0.242])
yticks(linspace(0.05, 0.25, 5))
xticks(bvals)
legend(Location="northwest")
ylabel('E-S signal contrast')
xlabel('b-value (s/mm^2)')
grid on
ax=gca();
ax.FontSize=12;

saveas(f2, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', 'ES_Contrast.png'))

% % Calculate effect sizes
% EffectSizes = (long_contrast-short_contrast)./sqrt(long_SE.^2+short_SE.^2);

%% Long/Short contrast

indx = 1;
S_long = signals(indx, 7:11, 1);
SE_long = signals(indx, 7:11, 2);
S_short = signals(indx, 2:6, 1);
SE_short = signals(indx, 2:6, 2);

diff = S_long-S_short;
diff_CI = 1.96*sqrt(SE_long.^2 + SE_short.^2);