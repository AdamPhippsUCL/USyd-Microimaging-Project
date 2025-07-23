% Display signal errors and modelling results in ROIs

clear;
projectfolder = pwd;


%% Sample, Image, and ROI

% Sample
SampleName = '20250224_UQ4'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

% ROI
ROIname = 'UQ4B Lesion 3+3'; 
% 'UQ4B Lesion 3+3'
% 'UQ4M Lesion 3+4'
% 'UQ4N Benign Glandular'
% 'UQ6N Lesion 4+4'
% 'UQ6M Stroma'

% Load ROI mask
ROI = load(fullfile(projectfolder, 'Scripts', 'ROIs', [ROIname '.mat'])).mask;
% ROI = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'BMASK.mat')).BMASK;

% Load composition and get ROI composition
COMPOSITION = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;
SampleMask = sum(COMPOSITION,4)>0;
ROI = ROI.*SampleMask;

ROI_COMP = reshape(COMPOSITION, [], 3);
flatMASK = ROI(:);
ROI_COMP = ROI_COMP(logical(flatMASK), :);

% Reshape for colouring markers
X = ROI_COMP;
X(:, [1,2]) = X(:, [2,1]);

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


% figure;
% bshift = 25;
% % Image
% for seriesindx = 2:length(SeriesDescriptions)

seriesindx = 11;
    
SeriesDescription = SeriesDescriptions{seriesindx};
bval = scheme(seriesindx).bval;

% Load normalized image
ImageArray = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, SeriesDescription, 'normalisedImageArray.mat')).ImageArray;

% Load signal measurements
signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'signals.mat')).signals;
signals = squeeze(signals(:,seriesindx,1));

% Load RESULTS
RESULTS = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'RESULTS.mat')).RESULTS;
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


figure
scatter( (pred_ROI+img_ROI)/2 , img_ROI-pred_ROI , 6, 'filled', 'MarkerFaceAlpha', 0.7, CData=X)
hold on
yline(LOA(1), '--')
yline(LOA(2))
yline(LOA(3))


%% Modelling results

ModelName = 'RDI - 2 compartment - 4 param (S0)';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, [ModelName ''], schemename, fittingtechnique);

fit_fIC = load(fullfile(outputfolder, 'fIC.mat')).fIC;
fit_fIC = fit_fIC(logical(ROI));
fit_dIC = load(fullfile(outputfolder, 'dIC.mat')).dIC;
fit_dIC = fit_dIC(logical(ROI));
fit_R = load(fullfile(outputfolder, 'R.mat')).R;
fit_R = fit_R(logical(ROI));
fit_dEES = load(fullfile(outputfolder, 'dEES.mat')).dEES;
fit_dEES = fit_dEES(logical(ROI));
fit_AIC = load(fullfile(outputfolder, 'AIC.mat')).AIC;
fit_AIC = fit_AIC(logical(ROI));

% Load ESL modelling estimates
ESL_Model_RESULTS = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'Modelling', 'RESULTS.mat')).RESULTS;
S_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'S'), strcmp({ESL_Model_RESULTS(:).ModelType}, ModelName))).ModelParams;
E_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'G'), strcmp({ESL_Model_RESULTS(:).ModelType}, ModelName))).ModelParams;
L_params = [1, 2]; % ADC parameters...


% Predicted model parameters
pred_fIC = S_params(1).*ROI_COMP(:,1) + E_params(1).*ROI_COMP(:,2);
pred_R = S_params(2).*ROI_COMP(:,1) + E_params(2).*ROI_COMP(:,2) + 6.5.*ROI_COMP(:,3);
pred_dIC = S_params(3).*ROI_COMP(:,1) + E_params(3).*ROI_COMP(:,2) + L_params(2).*ROI_COMP(:,3);
pred_dEES = S_params(4).*ROI_COMP(:,1) + E_params(4).*ROI_COMP(:,2) + L_params(2).*ROI_COMP(:,3);

% fIC
% m = max([pred_fIC; fit_fIC; 0.4]);
MAX = max([pred_fIC; fit_fIC]);
MIN = min([pred_fIC; fit_fIC]);

% figure
nexttile;
scatter(pred_fIC, fit_fIC, '*', CData=X)
hold on
plot([MIN,MAX],[MIN,MAX], '--', color = 'k')
xlim([MIN-0.02 MAX+0.02])
ylim([MIN-0.02 MAX+0.02])
xlabel('Predicted sphere fraction')
ylabel('Estimated sphere fraction')
grid on
xticks(linspace(0,1,11))
yticks(linspace(0,1,11))
title('Ball + Sphere Model')

% dEES
MAX = max([pred_dEES; fit_dEES]);
MIN = min([pred_dEES; fit_dEES]);
% figure
nexttile;
scatter(pred_dEES, fit_dEES, '*', CData=X)
hold on
plot([MIN,MAX],[MIN,MAX], '--', color = 'k')
xlim([MIN-0.05 MAX+0.05])
ylim([MIN-0.05 MAX+0.05])
grid on
xlabel('Predicted D_{out} (x10^{-3} mm^2/s)')
ylabel('Estimated D_{out} (x10^{-3} mm^2/s)')
grid on
xticks(linspace(0,3,16))
yticks(linspace(0,3,16))
title('Ball + Sphere Model')

saveas(f, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', ['ROI Results ' ROIname '.png']))

figure
scatter(pred_R, fit_R)