% MATLAB Script to average over 4th dimension components of an image 

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescription = '40u_DtiSE_2012_SPOIL10%';

% Use denoised data
UseDenoisedData = false;

%% Load image

switch UseDenoisedData
    case true
        % Load image array
        ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
    case false
        % Load image array
        ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
end


%% Aveerage over 4th dimension

avgImageArray = mean(ImageArray, 4);


% figure;
% imshow(squeeze(avgImageArray(60,:,:)), [])


%% Save 


switch UseDenoisedData
    case true
        save(fullfile(ImagingDataFolder, 'MAT DN' ,SampleName, SeriesDescription, 'avgImageArray.mat'), "avgImageArray");
    case false
        save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'avgImageArray.mat'), "avgImageArray");
end