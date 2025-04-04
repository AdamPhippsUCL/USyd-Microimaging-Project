% SCRIPT to auto generate tissue masks using 3DMGE, D, FA

clear;
projectfolder = pwd;

%% Sample and Image details

% Sample name
SampleName = '20250224_UQ4';

% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedDatas
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end

% 3D-MGE sequence
MGE_SeriesDescription = '3DMGE_20u';
MGE = load(fullfile(ImagingDataFolder, SampleName, MGE_SeriesDescription, 'avgImageArray.mat')).avgImageArray;

% dwFA
DTI_SeriesDescription = '40u_DtiSE_2012_SPOIL10% (20 micron)';
dwFA = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'dwFA.mat')).dwFA;


%% Generate masks

switch SampleName

    case '20250224_UQ4'
        
        STROMA = (dwFA>8e-5).*and(MGE<1.2e-7, MGE>5e-8);
        LUMEN = (MGE>1.2e-7);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>5e-8);

end


%% Display masks

displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = logical(GLANDULAR);
displaymasks(:,:,:,2) = logical(STROMA);
displaymasks(:,:,:,3) = logical(LUMEN);

figure
imshow(squeeze(MGE(128,:,:)),[]);
hold on
mask = imshow(squeeze(displaymasks(128,:,:,:)));
set(mask, 'AlphaData', 0.2)


%% Save masks

folder = fullfile(projectfolder, 'Outputs', 'Masks',SampleName, MGE_SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'LUMEN.mat'), 'LUMEN');
save(fullfile(folder, 'STROMA.mat'), 'STROMA');
save(fullfile(folder, 'GLANDULAR.mat'), 'GLANDULAR');