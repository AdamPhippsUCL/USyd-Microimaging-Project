% Script to downsample images

clear;
projectfolder = pwd;

%% Image details

% Imaging data folder 
ImagingDataFolder = fullfile(projectfolder, 'Imaging Data');

% Sample name
SampleName = '20250407_UQ5';

% Series description
SeriesDescriptions = {...
    % '3DMGE_20u'
    'SE_b0_SPOIL5%',...
    'STEAM_ShortDELTA_15',...
    'STEAM_ShortDELTA_20',...
    'STEAM_ShortDELTA_30',...
    'STEAM_ShortDELTA_40',...
    'STEAM_ShortDELTA_50',...
    'STEAM_LongDELTA_40',...
    'STEAM_LongDELTA_60',...
    'STEAM_LongDELTA_80',...
    'STEAM_LongDELTA_100',...
    'STEAM_LongDELTA_120'...
    };


% Downsample settings
downsamplewindow = [3,3,4]; 
overlap = false;

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    % Load image
    ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;

    % Load dinfo
    dinfo = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;

    % Downsample image
    ImageArray = downsample(ImageArray, downsamplewindow, overlap=overlap);

    % Update voxel size information in dinfo
    [dinfo(:).PixelSpacing] = deal([dinfo(1).PixelSpacing].*downsamplewindow(1:2));
    [dinfo(:).SliceThickness] = deal([dinfo(1).SliceThickness].*downsamplewindow(3));

    % Save downsampled data
    SeriesDescription = [SeriesDescription ' (DS)' ];
    folder = fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription);
    mkdir(folder)
    save(fullfile(folder, 'axialImageArray.mat'), 'ImageArray');
    save(fullfile(folder, 'axialdinfo.mat'), 'dinfo');

    try
        save(fullfile(folder, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');
    catch
        disp('')
    end

    Meta = struct();
    Meta.downsamplewindow = downsamplewindow;
    Meta.overlap = overlap;
    save(fullfile(folder, 'Meta.mat'), 'Meta');


end
