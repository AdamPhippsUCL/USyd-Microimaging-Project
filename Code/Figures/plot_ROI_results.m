% Display signal errors and modelling results in ROIs

clear;
projectfolder = pwd;


%% Sample, Image, and ROI

% Sample
SampleName = '20250224_UQ4'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

% ROI
ROIname = 'UQ4B and UQ4M Lesion 3+3'; 
% 'UQ4B Lesion 3+3'
% 'UQ4M Lesion 3+3'
% 'UQ4B and UQ4M Lesion 3+3'
% 'UQ6N Lesion 4+4'
% 'UQ6M Stroma'

% Load ROI mask
ROI = load(fullfile(projectfolder, 'Code', 'ROIs', [ROIname '.mat'])).mask;

% Load composition and get ROI composition
COMPOSITION = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;
SampleMask = sum(COMPOSITION,4)>0;
ROI = ROI.*SampleMask;

ROI_COMP = reshape(COMPOSITION, [], 3);
flatMASK = ROI(:);
ROI_COMP = ROI_COMP(logical(flatMASK), :);

X = ROI_COMP;
%% Signal error

% Scheme
scheme = load(fullfile(projectfolder, "Schemes", "20250224_UQ4 AllDELTA.mat")).scheme;

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


seriesindx = 11;
    
SeriesDescription = SeriesDescriptions{seriesindx};
bval = scheme(seriesindx).bval;
DELTA = scheme(seriesindx).DELTA;

% Load normalized image
ImageArray = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, SeriesDescription, 'normalisedImageArray.mat')).ImageArray;

% Load signal measurements
signals = load(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'signals.mat')).signals;
signals = squeeze(signals(:,seriesindx,1));

% Load RESULTS
RESULTS = load(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'RESULTS.mat')).RESULTS;
LOA = RESULTS(seriesindx).LOA;

% Predicted signals
signals = reshape(signals, [1,1,1,3]);
pred = sum(COMPOSITION.*repmat(signals, [size(COMPOSITION, 1:3)]), 4);

pred_ROI = pred(logical(SampleMask.*ROI));
img_ROI = ImageArray(logical(SampleMask.*ROI));

%     if seriesindx <=6
%         thisbshift = -bshift;
%     else
%         thisbshift = bshift;
%     end
%     scatter(bval*ones(1,numel(pred_ROI))+thisbshift, pred_ROI-img_ROI, '*', 'k')
%     hold on
% end
    
% m = max([pred_ROI; img_ROI; 0.6]);


% MAX = max([pred_ROI; img_ROI]);
% MIN = min([pred_ROI; img_ROI]);
% 
% f=figure;
% f.Position = [423   0   939   729];
% tiledlayout(2,2,'Padding', 'compact', 'TileSpacing', 'compact');
% nexttile;
% axis off
% 
% nexttile;
% scatter(pred_ROI, img_ROI, '*', CData=X)
% hold on
% plot([MIN,MAX],[MIN,MAX], '--', color = 'k')
% xlim([MIN-0.02 MAX+0.02])
% ylim([MIN-0.02 MAX+0.02])
% grid on
% xticks(linspace(0,1,11))
% yticks(linspace(0,1,11))
% xlabel('Predicted signal')
% ylabel('Measured signal')
% title('b = 2000 s/mm^2 , \Delta = 120 ms')



avg = (pred_ROI+img_ROI)/2;
diff = img_ROI-pred_ROI;

f1=figure;
scatter(avg, diff, 20, 'filled', 'MarkerFaceAlpha', 1, CData=X, HandleVisibility='off')
hold on
% yline(0, '--', HandleVisibility='off')
yline(LOA(1), '-', DisplayName = 'Bias (Benign)', LineWidth=1)
yline(LOA(2), '--', HandleVisibility='off', LineWidth=1)
yline(LOA(3), '--', DisplayName = '95% LOA (Benign)', LineWidth=1)
grid on
switch ROIname
    case 'UQ6N Lesion 4+4'
        ylim([-0.22 0.44])
        xlim([0.37, 0.68])
    case 'UQ4B Lesion 3+3'
        xlim([0.18, 0.58])
        ylim([-0.25 0.25])
    case 'UQ4B and UQ4M Lesion 3+3'
        xlim([0.18, 0.58])
        ylim([-0.25 0.25])
end
yticks(linspace(-1, 1, 21))
xticks(linspace(0, 1, 11))
% xlim([0.05*floor(min(avg)/0.05)-0.1 0.05*ceil(max(avg)/0.05)+0.1])
% ylim([-0.05*ceil(max(abs(diff))/0.05)-0.1 0.05*ceil(max(abs(diff))/0.05)+0.1])
legend(Location="northwest")

xlabel('Mean of Predicted and Measured Signal')
ylabel('Measured Signal - Predicted Signal')
title(['b = ' num2str(bval) 's/mm^2 ; Delta = ' num2str(DELTA) ' ms'])

ax1 = gca();
ax1.FontSize = 12;
f1.Position = [488   242   660   400];
saveas(f1, fullfile(projectfolder, 'Figures', ['Bland Altman ' ROIname ' b' num2str(bval) ' Delta' num2str(DELTA) '.png']))


%% Modelling results

ModelName = 'Ball+Sphere';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, [ModelName ''], schemename, fittingtechnique);

fit_fs = load(fullfile(outputfolder, 'fs.mat')).fs;
fit_fs = fit_fs(logical(ROI));
fit_Ds = load(fullfile(outputfolder, 'Ds.mat')).Ds;
fit_Ds = fit_Ds(logical(ROI));
fit_R = load(fullfile(outputfolder, 'R.mat')).R;
fit_R = fit_R(logical(ROI));
fit_Db = load(fullfile(outputfolder, 'Db.mat')).Db;
fit_Db = fit_Db(logical(ROI));
fit_AIC = load(fullfile(outputfolder, 'AIC.mat')).AIC;
fit_AIC = fit_AIC(logical(ROI));

% Load ESL modelling estimates
ESL_Model_RESULTS = load(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'Modelling', 'RESULTS.mat')).RESULTS;
S_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'S'), strcmp({ESL_Model_RESULTS(:).ModelName}, ModelName))).ModelParams;
E_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'E'), strcmp({ESL_Model_RESULTS(:).ModelName}, ModelName))).ModelParams;
L_params = [1, 2]; % ADC parameters...

% Predicted model parameters
pred_fs = S_params(1).*ROI_COMP(:,2) + E_params(1).*ROI_COMP(:,1);
pred_R = S_params(2).*ROI_COMP(:,2) + E_params(2).*ROI_COMP(:,1) + 6.5.*ROI_COMP(:,3);
pred_Ds = S_params(3).*ROI_COMP(:,2) + E_params(3).*ROI_COMP(:,1) + L_params(2).*ROI_COMP(:,3);
pred_Db = S_params(4).*ROI_COMP(:,2) + E_params(4).*ROI_COMP(:,1) + L_params(2).*ROI_COMP(:,3);

% Limits of agreement
load(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'Modelling', 'fs_LOA.mat'))
load(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'Modelling', 'Db_LOA.mat'))

% == BLAND ALTMAN

% Sphere fraction
fs_avg = (pred_fs+fit_fs)/2;
fs_diff = fit_fs-pred_fs;
f2=figure;
scatter(fs_avg, fs_diff, 20, 'filled', 'MarkerFaceAlpha', 1, CData=X, HandleVisibility='off')
hold on
yline(fs_LOA(1), '-', DisplayName = 'Bias (Benign)', LineWidth=1)
yline(fs_LOA(2), '--', HandleVisibility='off', LineWidth=1)
yline(fs_LOA(3), '--', DisplayName = '95% LOA (Benign)', LineWidth=1)
switch ROIname
    case 'UQ6N Lesion 4+4'
        xlim([0.23 0.42])
        ylim([-0.13, 0.32])
    case 'UQ4B Lesion 3+3'
        xlim([0.12 0.36])
        ylim([-0.16, 0.22])
    case 'UQ4B and UQ4M Lesion 3+3'
        xlim([0.12 0.36])
        ylim([-0.16, 0.22])
end
xticks(linspace(0, 1, 21))
yticks(linspace(-0.5, 0.5, 11))
xlabel('Mean of Predicted and Estimated Sphere Fraction')
ylabel('Estimated - Predicted Sphere Fraction')
legend(Location="northwest")
grid on
ax2 = gca();
ax2.FontSize = 12;
f2.Position = [488   242   660   400];


% Db 
Db_avg = (pred_Db+fit_Db)/2;
Db_diff = fit_Db-pred_Db;
f3=figure;
scatter(Db_avg, Db_diff, 20, 'filled', 'MarkerFaceAlpha', 1, CData=X, HandleVisibility='off')
hold on
yline(Db_LOA(1), '-', DisplayName = 'Bias (Benign)', LineWidth=1)
yline(Db_LOA(2), '--', HandleVisibility='off', LineWidth=1)
yline(Db_LOA(3), '--', DisplayName = '95% LOA (Benign)', LineWidth=1)
switch ROIname
    case 'UQ6N Lesion 4+4'
        xlim([0.33 0.72])
        ylim([-0.64, 0.64])
    case 'UQ4B Lesion 3+3'
        xlim([0.38 0.98])
        ylim([-0.52, 0.62])
    case 'UQ4B and UQ4M Lesion 3+3'
        xlim([0.36 1.04])
        ylim([-0.52, 0.62])
end
xticks(linspace(0, 3, 31))
yticks(linspace(-2, 2, 21))
xlabel('Mean of Predicted and Estimated D_{b} (x10^{-3} mm^2/s)')
ylabel('Estimated - Predicted D_{b} (x10^{-3} mm^2/s)')
legend(Location="northwest")
grid on
ax3 = gca();
ax3.FontSize = 12;
f3.Position = [488   242   660   400];

ax2.Position = ax3.Position;

saveas(f3, fullfile(projectfolder, 'Figures', ['Bland Altman_' ROIname ' Db.png']))
saveas(f2, fullfile(projectfolder, 'Figures', ['Bland Altman ' ROIname ' Sphere_Fraction.png']))

% 
% f=figure;
% ax = axes();
% ax.Position = ax1.Position;
% xticks([])
% yticks([])
% f.Position = [488   242   660   400];
% saveas(f, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', ['Histo_Positioner.png']))
% 



% fIC
% m = max([pred_fIC; fit_fIC; 0.4]);
% MAX = max([pred_fIC; fit_fIC]);
% MIN = min([pred_fIC; fit_fIC]);
% 
% figure
% nexttile;
% scatter(pred_fIC, fit_fIC, '*', CData=X)
% hold on
% plot([MIN,MAX],[MIN,MAX], '--', color = 'k')
% xlim([MIN-0.02 MAX+0.02])
% ylim([MIN-0.02 MAX+0.02])
% xlabel('Predicted sphere fraction')
% ylabel('Estimated sphere fraction')
% grid on
% xticks(linspace(0,1,11))
% yticks(linspace(0,1,11))
% title('Ball + Sphere Model')
% 
% dEES
% MAX = max([pred_dEES; fit_dEES]);
% MIN = min([pred_dEES; fit_dEES]);
% figure
% nexttile;
% scatter(pred_dEES, fit_dEES, '*', CData=X)
% hold on
% plot([MIN,MAX],[MIN,MAX], '--', color = 'k')
% xlim([MIN-0.05 MAX+0.05])
% ylim([MIN-0.05 MAX+0.05])
% grid on
% xlabel('Predicted D_{out} (x10^{-3} mm^2/s)')
% ylabel('Estimated D_{out} (x10^{-3} mm^2/s)')
% grid on
% xticks(linspace(0,3,16))
% yticks(linspace(0,3,16))
% title('Ball + Sphere Model')
% 
% saveas(f, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', ['ROI Results ' ROIname '.png']))
