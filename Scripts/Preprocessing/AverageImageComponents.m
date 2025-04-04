% Script to average over 4th dimension component of an image 

clear;
projectfolder = pwd;

%% Image details

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescription =  'SE_b0_SPOIL5%';

% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end

%% Load image

% Load image array
ImageArray = load(fullfile(ImagingDataFolder, SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;


%% Aveerage over 4th dimension

avgImageArray = mean(ImageArray, 4);

%% Save 

save(fullfile(ImagingDataFolder, SampleName, SeriesDescription, 'avgImageArray.mat'), "avgImageArray");
