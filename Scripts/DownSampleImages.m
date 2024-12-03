% MATLAB script to downsample images
clear;

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


% Downsample settings
downsamplewindow = [2,2,2]; 
overlap = false;

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    % Load image
    ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;

    % Load dinfo
    dinfo = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;

    % Load VoxelCoordinates
    VoxelCoordinates = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialVoxelCoordinates.mat')).VoxelCoordinates;
 

    % Downsample image
    ImageArray = downsample(ImageArray, downsamplewindow, overlap=overlap);

    % Downsample voxel coordinates
    VoxelCoordinates = downsample(VoxelCoordinates, downsamplewindow, overlap=overlap);

    % Update voxel size information in dinfo
    [dinfo(:).PixelSpacing] = deal([dinfo(1).PixelSpacing].*downsamplewindow(1:2));
    [dinfo(:).SliceThickness] = deal([dinfo(1).SliceThickness].*downsamplewindow(3));


    % Save downsampled data
    SeriesDescription = [SeriesDescription ' (' num2str(dinfo(1).SliceThickness*1000) ' micron)' ];
    folder = fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription);
    mkdir(folder)
    save(fullfile(folder, 'axialImageArray.mat'), 'ImageArray');
    save(fullfile(folder, 'axialdinfo.mat'), 'dinfo');
    save(fullfile(folder, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');

    Meta = struct();
    Meta.downsamplewindow = downsamplewindow;
    Meta.overlap = overlap;
    save(fullfile(folder, 'Meta.mat'), 'Meta');


end
