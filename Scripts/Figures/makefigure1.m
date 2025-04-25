% Script to make figure 1

clear;
projectfolder = pwd;


% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end


MGE_SeriesDescription = '3DMGE_20u';
DTI_SeriesDescription = '40u_DtiSE_2012_SPOIL10% (20 micron)';


%% Left column details

SampleName = '20250224_UQ4';

% Load images
MGE1 = load(fullfile(ImagingDataFolder, SampleName, MGE_SeriesDescription, 'avgImageArray.mat')).avgImageArray;
dwFA1 = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'dwFA.mat')).dwFA;

% Load masks
maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', SampleName, MGE_SeriesDescription);
GLANDULAR1 = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA1 = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN1 = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;



%% Left column

displaymasks1 = zeros([size(GLANDULAR1), 3]);
displaymasks1(:,:,:,1) = logical(GLANDULAR1);
displaymasks1(:,:,:,2) = logical(STROMA1);
displaymasks1(:,:,:,3) = logical(LUMEN1);

sl=140;
cols = 400:630;
rows = 54:200;

f1=figure;
imshow(squeeze(MGE1(sl,rows,cols)),[0 prctile(squeeze(MGE1(sl,rows,cols)), 99.9, 'all')]);
cb=colorbar;
cb.Label.String = 'MGE signal intensity';
title('MGE')
ax = gca();
ax.FontSize = 10;
ax.Position = [0.06         0.18         0.79         0.68];

f2=figure;
imshow(squeeze(dwFA1(sl,rows,cols)),[0 prctile(squeeze(dwFA1(sl,rows,cols)), 99.9, 'all')]);
cb=colorbar;
cb.Label.String = 'D x FA (mm/s^2)';
title('D x FA')
ax = gca();
ax.FontSize = 10;
ax.Position = [0.06         0.18         0.79         0.68];


f3=figure;
imshow(squeeze(MGE1(sl,rows,cols)),[0 prctile(squeeze(MGE1(sl,rows,cols)), 99.9, 'all')]);
hold on
mask = imshow(squeeze(displaymasks1(sl,rows,cols,:)));
set(mask, 'AlphaData', 0.2)
title('ESL segmentation')
ax = gca();
ax.FontSize = 10;
ax.Position = [0.06         0.18         0.79         0.68];

saveas(f1, fullfile(projectfolder, "Scripts", "Figures", "MGE1.fig"))
saveas(f1, fullfile(projectfolder, "Scripts", "Figures", "MGE1.png"))
saveas(f2, fullfile(projectfolder, "Scripts", "Figures", "DFA1.fig"))
saveas(f2, fullfile(projectfolder, "Scripts", "Figures", "DFA1.png"))
saveas(f3, fullfile(projectfolder, "Scripts", "Figures", "ESL1.fig"))
saveas(f3, fullfile(projectfolder, "Scripts", "Figures", "ESL1.png"))

%% Right column details

SampleName = '20250224_UQ4';

% Load images
MGE2 = load(fullfile(ImagingDataFolder, SampleName, MGE_SeriesDescription, 'avgImageArray.mat')).avgImageArray;
dwFA2 = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'dwFA.mat')).dwFA;

% Load masks
maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', SampleName, MGE_SeriesDescription);
GLANDULAR2 = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA2 = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN2 = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;


%% Right column

displaymasks2 = zeros([size(GLANDULAR2), 3]);
displaymasks2(:,:,:,1) = logical(GLANDULAR2);
displaymasks2(:,:,:,2) = logical(STROMA2);
displaymasks2(:,:,:,3) = logical(LUMEN2);

sl=115;
cols = 200:400;
rows = 55:182;

f1=figure;
imshow(squeeze(MGE2(sl,rows,cols)),[0 prctile(squeeze(MGE2(sl,rows,cols)), 99.9, 'all')]);
cb=colorbar;
cb.Label.String = 'MGE signal intensity';
title('MGE')
ax = gca();
ax.FontSize = 10;
ax.Position = [0.06         0.18         0.79         0.68];

f2=figure;
imshow(squeeze(dwFA2(sl,rows,cols)),[0 prctile(squeeze(dwFA2(sl,rows,cols)), 99.9, 'all')]);
cb=colorbar;
cb.Label.String = 'D x FA (mm/s^2)';
title('D x FA')
ax = gca();
ax.FontSize = 10;
ax.Position = [0.06         0.18         0.79         0.68];


f3=figure;
imshow(squeeze(MGE2(sl,rows,cols)),[0 prctile(squeeze(MGE2(sl,rows,cols)), 99.9, 'all')]);
hold on
mask = imshow(squeeze(displaymasks2(sl,rows,cols,:)));
set(mask, 'AlphaData', 0.2)
title('ESL segmentation')
ax = gca();
ax.FontSize = 10;
ax.Position = [0.06         0.18         0.79         0.68];


tissue_labels = {'S', 'G', 'L'};

saveas(f1, fullfile(projectfolder, "Scripts", "Figures", "MGE2.fig"))
saveas(f1, fullfile(projectfolder, "Scripts", "Figures", "MGE2.png"))
saveas(f2, fullfile(projectfolder, "Scripts", "Figures", "DFA2.fig"))
saveas(f2, fullfile(projectfolder, "Scripts", "Figures", "DFA2.png"))
saveas(f3, fullfile(projectfolder, "Scripts", "Figures", "ESL2.fig"))
saveas(f3, fullfile(projectfolder, "Scripts", "Figures", "ESL2.png"))