% MATLAB SCRIPT TO RESIZE DTI IMAGE

clear;

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescription = '40u_DtiSE_2012_SPOIL10%';

% Use denoised data
UseDenoisedData = true;

%% Load image

switch UseDenoisedData
    case true
        % Load image array
        ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
        % Load dinfo
        dinfo=load(fullfile(ImagingDataFolder, 'MAT DN' ,SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;
    case false
        % Load image array
        ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
        % Load dinfo
        dinfo=load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;
end




%% Resize

ResFactor = [2,2,2];
newImageArray = zeros([ResFactor.*size(ImageArray,1:3) size(ImageArray,4)]);

% Resize and apply Gaussian smoothing
for t = 1:size(ImageArray, 4)
    newdata = imresize3(ImageArray(:,:,:,t), ResFactor.*size(ImageArray,1:3), 'linear');
    newImageArray(:,:,:,t) = imgaussfilt(newdata, 0.5);
    % newImageArray(:,:,:,t) = newdata;
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


% Make directory
switch UseDenoisedData
    case true
        folder = fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription);
        mkdir(folder)
    case false
        folder = fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription);
        mkdir(folder)
end

% save(fullfile(folder, 'axialImageArray.mat'), 'ImageArray');
save(fullfile(folder, 'axialImageArray.mat'), 'ImageArray', '-v7.3');
save(fullfile(folder, 'axialdinfo.mat'), 'dinfo');