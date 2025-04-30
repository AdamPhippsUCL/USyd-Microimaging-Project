% Script to make figure 2

clear;
projectfolder = pwd;


outputfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement'); 

samplename = 'ESMRMB';
signals = load(fullfile(outputfolder, samplename, 'signals.mat')).signals;
scheme = load(fullfile(outputfolder, samplename, 'scheme.mat')).scheme;
nscheme = length(scheme);


%% Figure


f=figure;
bshift = 12;

% Short Delta

indices = 2:6;
s = signals(:,indices,:);
bvals = [scheme(indices).bval];

errorbar(bvals-1*bshift, s(2,:,1), s(2,:,3), '-*', color=[0.8500 0.3250 0.0980], DisplayName='E (Short \Delta)');
hold on
errorbar(bvals-1*bshift, s(1,:,1), s(1,:,3), '-*', color=[0.4660 0.6740 0.1880], DisplayName='S (Short \Delta)');
errorbar(bvals-1*bshift, s(3,:,1), s(3,:,3), '-*', color=[0 0.4470 0.7410], DisplayName = 'L (Short \Delta)');

% Long Delta
indices = 7:11;
s = signals(:,indices,:);

errorbar(bvals+1*bshift, s(2,:,1), s(2,:,3), '-.*', color=[0.8500 0.3250 0.0980], DisplayName='E (Long \Delta)');
errorbar(bvals+1*bshift, s(1,:,1), s(1,:,3), '-.*', color=[0.4660 0.6740 0.1880], DisplayName='S (Long \Delta)');
errorbar(bvals+1*bshift, s(3,:,1), s(3,:,3), '-.*', color=[0 0.4470 0.7410], DisplayName = 'L (Long \Delta)');

xticks(bvals); 
xticklabels(bvals)
ylim([-0.0,0.7])
ylabel('Normalised signal estimate')
xlim([800,2200])
xlabel('b-value (s/mm^{2})')
xticks(bvals)
xticklabels(["1000", "1250", "1500", "1750", "2000"])
grid on
legend('NumColumns', 2);
ax=gca();
ax.FontSize=12;
f.Position = [488   242   700   480];

% Save figure
saveas(f, fullfile(projectfolder, "Scripts", "ESMRMB", "Figures", "Figure2.fig"))
saveas(f, fullfile(projectfolder, "Scripts", "ESMRMB", "Figures", "Figure2.png"))