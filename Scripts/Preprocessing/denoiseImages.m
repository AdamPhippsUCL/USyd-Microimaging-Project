% MATLAB Script to apply MP-PCA denoising to images

clear;
projectfolder = pwd;

%% Image details

% Imaging data folder 
ImagingDataFolder = fullfile(projectfolder, 'Imaging Data');

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescriptions = {...
     'SE_b0_SPOIL5%',...
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



