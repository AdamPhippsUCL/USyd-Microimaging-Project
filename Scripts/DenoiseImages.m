% MATLAB Script to apply MP-PCA denoising to images

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescriptions = {...
    % 'SE_b0_SPOIL10%',...
    'STEAM_ShortDELTA_15',...
    'STEAM_ShortDELTA_20',...
    'STEAM_ShortDELTA_30',...
    'STEAM_ShortDELTA_40',...
    'STEAM_ShortDELTA_50',...
    };


for seriesindx = 1:length(SeriesDescriptions)

    
    SeriesDescription = SeriesDescriptions{seriesindx}
    
    % Load image array
    ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
    
    % Load dinfo
    dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;
    
    try
        % Load voxel coordinates
        VoxelCoordinates = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'axialVoxelCoordinates.mat')).VoxelCoordinates;
    catch
        disp('')
    end

    % SELECT FRAMES FOR DENOISING
    ImageArray = ImageArray;
    dinfo = dinfo;

    %% Denoising
    
    % Define window
    window = [5 5 1];
    
    % Apply denoising
    ImageArrayDN = denoise(ImageArray, window);
    
    
    
    % Example slices
    % slice = 15;
    % figure;
    % imshow(ImageArray(:,:,slice,6), [0 5e-4])
    % figure;
    % imshow(ImageArrayDN(:,:,slice,6), [0 5e-4])
    
    
    %% Save results
    
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



