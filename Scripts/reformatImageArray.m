% MATLAB script to reformat shape of image arrays from Bruker system

%% Specify image

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241213_WSU1';

% Series description
SeriesDescription = '100001';


%% Load image array and DICOM 

ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'ImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'dinfo.mat')).dinfo;
VoxelCoordinates = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'VoxelCoordinates.mat')).VoxelCoordinates;

H = dinfo(1).Height;
W = dinfo(1).Width;


%% Separate acquisitions

% OLD (works for UQ)

% slices = unique([dinfo(:).sl]);
% AcquisitionNumbers = unique([dinfo(:).AcquisitionNumber]);
% 
% newImageArray = zeros([H, W, length(slices), length(AcquisitionNumbers)]);
% 
% acqfirstindices = [];
% 
% for acqindx = 1:length(AcquisitionNumbers)
% 
%     acqnum = AcquisitionNumbers(acqindx);
% 
%     framebools = [dinfo(:).AcquisitionNumber] == acqnum;
%     frameindices = find(framebools);
%     acqfirstindices = [acqfirstindices, frameindices(1)];
% 
%     thisimg = ImageArray(:,:,framebools);
%     newImageArray(:,:,:,acqnum) = thisimg;
% 
% end
% 
% newdinfo = dinfo(acqfirstindices);
% 
% % remove slice dependent dinfo fields
% dinfo = rmfield(newdinfo, {'TinSeries', 'sl', 'ImagePositionPatient', 'ImageOrientationPatient', 'Height', 'Width'});



% NEW

slices = unique([dinfo(:).sl]);

newImageArray = zeros([H, W, length(slices), length(dinfo)/length(slices)]);

for slindx = 1:length(slices)

    slice = slices(slindx);

    slicebools = [dinfo(:).sl] == slice;
    sliceindices = find(slicebools);

    thisimg = ImageArray(:,:,slicebools);
    newImageArray(:,:,slindx,:) = thisimg;

    if slindx == 1
        newdinfo = dinfo(sliceindices);
    end

end

% remove slice dependent dinfo fields
dinfo = rmfield(newdinfo, {'TinSeries', 'sl', 'ImagePositionPatient', 'ImageOrientationPatient', 'Height', 'Width'});



%% Change slice orientation to axial

if size(newImageArray,1)~=size(newImageArray,2)
    ImageArray = permute(newImageArray, [2,3,1,4]);
    VoxelCoordinates = permute(VoxelCoordinates, [2,3,1,4]);
else
    ImageArray = newImageArray;
end


%% Save reformated

save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat'), 'ImageArray');
save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');
save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialdinfo.mat'), 'dinfo');








