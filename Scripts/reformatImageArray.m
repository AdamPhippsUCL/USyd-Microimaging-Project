% MATLAB script to reformat shape of image arrays from Bruker system

%% Specify image

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241128_UQ1';

% Series description
SeriesDescription = '40u_verdict_seq5_v2';

%% Load image array and DICOM 

ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'ImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'dinfo.mat')).dinfo;
VoxelCoordinates = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'VoxelCoordinates.mat')).VoxelCoordinates;

H = dinfo(1).Height;
W = dinfo(1).Width;


%% Separate acquisitions

slices = unique([dinfo(:).sl]);
AcquisitionNumbers = unique([dinfo(:).AcquisitionNumber]);

newImageArray = zeros([H, W, length(slices), length(AcquisitionNumbers)]);

acqfirstindices = [];

for acqindx = 1:length(AcquisitionNumbers)
    
    acqnum = AcquisitionNumbers(acqindx);

    framebools = [dinfo(:).AcquisitionNumber] == acqnum;
    frameindices = find(framebools);
    acqfirstindices = [acqfirstindices, frameindices(1)];

    thisimg = ImageArray(:,:,framebools);
    newImageArray(:,:,:,acqnum) = thisimg;

end

newdinfo = dinfo(acqfirstindices);

% remove slice dependent dinfo fields
dinfo = rmfield(newdinfo, {'TinSeries', 'sl', 'ImagePositionPatient', 'ImageOrientationPatient', 'Height', 'Width'});


%% Change slice orientation to axial

ImageArray = permute(newImageArray, [2,3,1,4]);
VoxelCoordinates = permute(VoxelCoordinates, [2,3,1,4]);


%% Save reformated

save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat'), 'ImageArray');
save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');
save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialdinfo.mat'), 'dinfo');








