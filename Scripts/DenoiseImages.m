% MATLAB Script to apply MP-PCA denoising to images

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241128_UQ1';

% Series description
SeriesDescriptions = {...
    '40u_verdict_seq1_v2',...
    '40u_verdict_seq2_v2',...
    '40u_verdict_seq3_v2',...
    '40u_verdict_seq4_v2',...
    '40u_verdict_seq5_v2',...
    };


for seriesindx = 1:length(SeriesDescriptions)

    
    SeriesDescription = SeriesDescriptions{seriesindx}
    
    % Load image array
    ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
    
    % Load dinfo
    dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;
    
    % Load voxel coordinates
    VoxelCoordinates = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'axialVoxelCoordinates.mat')).VoxelCoordinates;
    
    %% Denoising
    
    % Define window
    window = [5 5];
    
    % Apply denoising
    ImageArrayDN = denoise(ImageArray, window);
    
    
    
    % % Example slices
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
    save(fullfile(folder, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');
    save(fullfile(folder, 'window.mat'), 'window');
    


end



