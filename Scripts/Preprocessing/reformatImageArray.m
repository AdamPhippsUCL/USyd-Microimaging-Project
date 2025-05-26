% Script to reformat shape of image arrays from Bruker system

clear;
projectfolder = pwd;

%% Specify sample and image

% Imaging data folder 
ImagingDataFolder = fullfile(projectfolder, 'Imaging Data');

% Sample name
SampleName = '20250524_UQ9';

% Series description
SeriesDescription =  'STEAM_DELTA_50';


%% Load image array and DICOM 

ImageArray = load(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'ImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'dinfo.mat')).dinfo;

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



%% Change slice orientation to axial

if size(newImageArray,1)~=size(newImageArray,2)
    ImageArray = permute(newImageArray, [2,3,1,4]);

else
    ImageArray = newImageArray;
end


%% Save reformated image and dinfo

save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialImageArray.mat'), 'ImageArray');
save(fullfile(ImagingDataFolder, 'MAT' ,SampleName, SeriesDescription, 'axialdinfo.mat'), 'dinfo');





