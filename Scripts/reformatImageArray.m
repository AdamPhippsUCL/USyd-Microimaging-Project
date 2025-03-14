% MATLAB script to reformat shape of image arrays from Bruker system

clear;

%% Specify image

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescription = '3DMGE_20u';


%% Load image array and DICOM 

ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'ImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'dinfo.mat')).dinfo;

try
    VoxelCoordinates = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'VoxelCoordinates.mat')).VoxelCoordinates;
catch
    disp('')
end

H = dinfo(1).Height;
W = dinfo(1).Width;

%% Separate acquisitions

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


% %% Resize array for CPC 20250131 data
% 
% oldsize = size(newImageArray);
% 
% temparray = zeros([64,64,7,4]);
% % Resize stack-wise
% for stindx = 1:size(newImageArray, 4)
%     temparray(:,:,:,stindx) = imresize3( ...
%         squeeze(newImageArray(:,:,:,stindx)),...
%         [64, 64, 7]...
%     );
% 
% end
% 
% newImageArray = temparray;


%% Change slice orientation to axial

if size(newImageArray,1)~=size(newImageArray,2)
    ImageArray = permute(newImageArray, [2,3,1,4]);

    try
        VoxelCoordinates = permute(VoxelCoordinates, [2,3,1,4]);
    catch
        disp('')
    end

else
    ImageArray = newImageArray;
end


%% Save reformated

save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat'), 'ImageArray');
save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialdinfo.mat'), 'dinfo');

try
    save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialVoxelCoordinates.mat'), 'VoxelCoordinates');
catch
    disp('')
end






