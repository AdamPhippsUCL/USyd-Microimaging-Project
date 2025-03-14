% MATLAB Script to save DICOM image as MAT array with corresponding dinfo structure

clear;

%% Initial definitions

% Sample name
SampleName = '20250224_UQ4';

% DICOM folder
DICOMfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\Imaging Data\20250224_101912_RB_Q4_RB_Q4_1_1\48\pdata\1\dicom";

% Imaging data folder (to save MAT images)
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Image type
imgtype = 'DTI';

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

if strcmp(SeriesDescription, '')
    SeriesDescription = num2str(dinfo(1).SeriesNumber);
end

% % USING NUMBER!!
% SeriesDescription = num2str(dinfo(1).SeriesNumber);

% PV Scaling
RS = [dinfo.RescaleSlope];
RI = [dinfo.RescaleIntercept];
SS = [dinfo.Private_2005_100e];


%% Diffusion Information (read manually from Method file)

switch imgtype 

    case 'DTI'

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

folder = fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription);
mkdir(folder);
save(fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription, 'dinfo.mat'), 'dinfo');
save(fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription, 'ImageArray.mat'), 'ImageArray');


% %% Voxel coordinates
% VoxelCoordinates = constructVoxelCoordinates(dinfo);
% save(fullfile(ImagingDataFolder, 'MAT' , SampleName, SeriesDescription, 'VoxelCoordinates.mat'), 'VoxelCoordinates');

