% Script to downsample images

clear;
projectfolder = pwd;

%% Image details

% Imaging data folder 
ImagingDataFolder = fullfile(projectfolder, 'Imaging Data');

% Sample name
SampleName = '20250524_UQ9'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

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


% Downsample settings (old function)
downsamplewindow = [3,3,2]; 
overlap = false;

% % Using imresize3
% imsize = [10,10,20];

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    % Load image
    ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;

    % Load dinfo
    dinfo = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;

    % Downsample image
    ImageArray = downsample(ImageArray, downsamplewindow, overlap=overlap);
    
    % % Using imresize3
    % n = size(ImageArray,4);
    % ImageArray_DS = zeros([imsize, n]);
    % for t = 1:n
    %     ImageArray_DS(:,:,:,t) = imresize3(ImageArray(:,:,:,t), [10,10,20], "AntiAliasing",true);
    % end
    % ImageArray = ImageArray_DS;

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
