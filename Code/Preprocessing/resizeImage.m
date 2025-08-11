% Script to resize DTI images

clear;
projectfolder = pwd;

%% Image details

% Sample name
SampleName = '20250524_UQ9';

% Series description
SeriesDescription = '40u_DtiSE_2012_SPOIL10%';

% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end

%% Load image and dinfo

% Load image array
ImageArray = load(fullfile(ImagingDataFolder, SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
% Load dinfo
dinfo=load(fullfile(ImagingDataFolder, SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;


%% Resize

ResFactor = [2,2,2];
newImageArray = zeros([ResFactor.*size(ImageArray,1:3) size(ImageArray,4)]);

% Resize and apply Gaussian smoothing
for t = 1:size(ImageArray, 4)
    newdata = imresize3(ImageArray(:,:,:,t), ResFactor.*size(ImageArray,1:3), 'linear');
    newImageArray(:,:,:,t) = imgaussfilt(newdata, 0.5);
end


%% Save new image

oldImageArray = ImageArray;
ImageArray = newImageArray;
clear newImageArray;

% Update voxel sizes
[dinfo(:).PixelSpacing] = deal([dinfo(1).PixelSpacing]./ResFactor(1:2));
[dinfo(:).SliceThickness] = deal([dinfo(1).SliceThickness]./ResFactor(3));


% Save downsampled data
SeriesDescription = [SeriesDescription ' (' num2str(dinfo(1).SliceThickness*1000) ' micron)' ];

folder = fullfile(ImagingDataFolder, SampleName, SeriesDescription);
mkdir(folder)

save(fullfile(folder, 'axialImageArray.mat'), 'ImageArray', '-v7.3');
save(fullfile(folder, 'axialdinfo.mat'), 'dinfo');