% MATLAB test script to generate mask using 3DMGE sequence

clear;

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescription = '3DMGE_20u';

% Use denoised data
UseDenoisedData = false;

%% Load image

switch UseDenoisedData
    case true
        % Load image array
        ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN' , SampleName, SeriesDescription, 'avgImageArray.mat')).avgImageArray;
    case false
        % Load image array
        ImageArray = load(fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription, 'avgImageArray.mat')).avgImageArray;
end


%% Testing on single slice

slice = squeeze(ImageArray(120,:,:));

%% Define thresholds

thres1 = 8.5e-8;
thres2 = 1.2e-7;

% STROMA
STROMA = (slice<thres1);
EPITHELIUM = and(slice>thres1, slice<thres2);
LUMEN = slice>thres2;


%% Generate ESL maps

ResFactor = 8;

MAP = zeros((1/ResFactor)*size(slice));
szmap = size(MAP);

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
    
        E = sum(EPITHELIUM( ...
            (rindx-1)*ResFactor+1 : rindx*ResFactor,...
            (cindx-1)*ResFactor+1 : cindx*ResFactor...
                ),...
                "all"...
        );

        MAP(rindx, cindx) = E./(ResFactor^2);

    end
end