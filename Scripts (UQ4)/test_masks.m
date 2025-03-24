% MATLAB SCRIPT to auto generate tissue masks using 3DMGE, D, FA

clear;

% Folders
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

% Sample name
SampleName = '20250224_UQ4';

% == 3DMGE sequence

MGE_SeriesDescription = '3DMGE_20u';
MGE = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, MGE_SeriesDescription, 'avgImageArray.mat')).avgImageArray;


% == FA and D

% Series description
DTI_SeriesDescription = '40u_DtiSE_2012_SPOIL10% (20 micron)';
% D = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'D.mat')).D;
% FA = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'FA.mat')).FA;
dwFA = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'dwFA.mat')).dwFA;


%% Generate masks

% STROMA
% STROMA = or(FA>0.15,D>0.015).*and(MGE<1.2e-7, MGE>4e-8);
% STROMA = (FA>0.15).*and(MGE<1.2e-7, MGE>4e-8);
STROMA = (dwFA>0.0016).*and(MGE<1.2e-7, MGE>5e-8);

% LUMEN
LUMEN = (MGE>1.2e-7);

% GLANDULAR
GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>5e-8);


%% Display masks

displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = logical(GLANDULAR);
displaymasks(:,:,:,2) = logical(STROMA);
displaymasks(:,:,:,3) = logical(LUMEN);

figure
imshow(squeeze(MGE(163,:,:)),[]);
hold on
mask = imshow(squeeze(displaymasks(163,:,:,:)));
set(mask, 'AlphaData', 0.2)


%% Save masks

folder = fullfile(projectfolder, 'Outputs', 'Masks',SampleName, MGE_SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'LUMEN.mat'), 'LUMEN');
save(fullfile(folder, 'STROMA.mat'), 'STROMA');
save(fullfile(folder, 'GLANDULAR.mat'), 'GLANDULAR');