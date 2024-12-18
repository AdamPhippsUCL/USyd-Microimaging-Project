% Test script to load DICOM and produce MAT image array with structure
% outling details of each frame.

% DICOM folder
DICOMfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\Imaging Data\20241128_112522_RB_Q1_RB_Q1_1_1\36\pdata\1\dicom";

% Imaging Data folder
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241128_UQ1';

%% Read DICOM information

% All DICOM filenames in folder
x = dir(DICOMfolder);
fnames = cell2mat(transpose({char(x(3:end).name)}));
dfolders = repmat([char(DICOMfolder) '/' ], [size(fnames,1), 1]);
dfnames = [string([dfolders  fnames])];
dfnames = cellstr(transpose(dfnames));

dinfo=dfparse(dfnames);

% Series description
SeriesDescription = dinfo(1).SeriesDescription;


% PV Scaling
RS = [dinfo.RescaleSlope];
RI = [dinfo.RescaleIntercept];
SS = [dinfo.Private_2005_100e];



%% Diffusion information

% Load DTI information structure
DTIstruct = load( fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription, 'DTIstruct.mat') ).DTIstruct;
bvals = DTIstruct.DiffusionBValue;
effbvals = DTIstruct.DiffusionEffBValue;

[dinfo(:).DiffusionBValue] = deal(0);
[dinfo(:).DiffusionEffBValue] = deal(0);

% Append to dinfo structure
for sl = unique([dinfo(:).sl])
    framebools = ([dinfo(:).sl] == sl);
    indices = find(framebools);
    for bindx = 1:length(indices)
        bval = bvals(bindx);
        indice = indices(bindx);
        dinfo(indice).DiffusionBValue = bvals(bindx);
        dinfo(indice).DiffusionEffBValue = effbvals(bindx);
    end
end


%% Construct Image Array

Nframes = length(dinfo);
ImageArray = zeros([dinfo(1).Height, dinfo(1).Width, Nframes]);

for frameindx = 1:Nframes
    img = double(dicomread(dfnames{frameindx}));
    img = (img*RS(frameindx)+RI(frameindx))/SS(frameindx);
    ImageArray(:,:,frameindx)=img;
end


%% Save dinfo structure and image array

save(fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription, 'dinfo.mat'), 'dinfo');
save(fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription, 'ImageArray.mat'), 'ImageArray');










% % Read from method file for diffusion b value and directions..
% fmethod = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\Imaging Data\20241128_112522_RB_Q1_RB_Q1_1_1\36\method";
% 
% fileread(fmethod)
% 
% 
% % 2dseq file
% f2dseq = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\Imaging Data\20241128_112522_RB_Q1_RB_Q1_1_1\36\pdata\1\2dseq";
% f1 = fopen(f2dseq,'r');
% A=fread(f1,'int16');


% dicomread(dfnames)