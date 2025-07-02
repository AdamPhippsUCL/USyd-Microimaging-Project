% MATLAB script to show SGL masks

clear;
projectfolder = pwd;

% Sample
SampleNum = 5;
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'};
SampleName = SampleNames{SampleNum};

% MGE image
MGE_seriesdescription = '3DMGE_20u';
MGE = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, MGE_seriesdescription, 'avgImageArray.mat')).avgImageArray;

% dwFA
dwFA = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', '40u_DtiSE_2012_SPOIL10% (20 micron)', 'dwFA.mat')).dwFA;

% Mask folder
maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', SampleName, MGE_seriesdescription);

% LOAD SGL MASKs
STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;

% Create 4D mask (color coded)
displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = logical(GLANDULAR);
displaymasks(:,:,:,2) = logical(STROMA);
displaymasks(:,:,:,3) = logical(LUMEN);


%% Define slice and region to present

xs = 1:240;
ys = 120;
zs = 1:640;

figure
imshow(squeeze(MGE(ys,xs,zs)),[0 prctile(squeeze(MGE(ys,xs,zs)), 99.9, 'all')]);
colorbar;

figure
imshow(squeeze(dwFA(ys,xs,zs)),[0 prctile(squeeze(dwFA(ys,xs,zs)), 99.9, 'all')]);

figure
imshow(squeeze(MGE(ys,xs,zs)),[0 prctile(squeeze(MGE(ys,xs,zs)), 99.9, 'all')]);
hold on
mask = imshow(squeeze(displaymasks(ys,xs,zs,:)));
set(mask, 'AlphaData', 0.2)