% MATLAB Script to apply MP-PCA denoising to images

clear;
projectfolder = pwd;

%% Image details

% Imaging data folder 
ImagingDataFolder = fullfile(projectfolder, 'Imaging Data');

% Sample name
SampleName = '20250524_UQ9';

% Series description
SeriesDescriptions = {...
     '40u_DtiSE_2012_SPOIL10%',...
     % 'SE_b0_SPOIL5%',...
     % 'STEAM_ShortDELTA_15',...
     % 'STEAM_ShortDELTA_20',...
     % 'STEAM_ShortDELTA_30',...
     % 'STEAM_ShortDELTA_40',...    
     % 'STEAM_ShortDELTA_50',...
     % 'STEAM_LongDELTA_40',...
     % 'STEAM_LongDELTA_60',...
     % 'STEAM_LongDELTA_80',...
     % 'STEAM_LongDELTA_100',...
     % 'STEAM_LongDELTA_120',...
     };


%% Denoising

% Define window
window = [5 5 1];

for seriesindx = 1:length(SeriesDescriptions)

    
    SeriesDescription = SeriesDescriptions{seriesindx}
    
    % Load image array
    ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
    
    % Load dinfo
    dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;

    % % SELECT FRAMES FOR DENOISING
    % ImageArray = ImageArray;
    % dinfo = dinfo;
    
    % Apply denoising
    ImageArrayDN = denoise(ImageArray, window);
        
    % Save results
    ImageArray = ImageArrayDN;
    
    folder = fullfile(ImagingDataFolder, 'MAT DN' ,SampleName, SeriesDescription);
    mkdir(folder);
    
    save(fullfile(folder, 'axialImageArray.mat'), 'ImageArray');
    save(fullfile(folder, 'axialdinfo.mat'), 'dinfo');
    try
        save(fullfile(folder, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');
    catch
        disp('')
    end    
    save(fullfile(folder, 'window.mat'), 'window');
    

end



