% Script to define ROI over low resolution images

clear;
projectfolder = pwd;

%% Sample and image

% Sample
SampleName = '20250224_UQ4'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

% Image
SeriesDescription = 'STEAM_LongDELTA_120 (DS)';

% Load normalized image
ImageArray = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, SeriesDescription, 'normalisedImageArray.mat')).ImageArray;

% Base image
BaseSeriesDescription = '3DMGE_20u';

%

% Load normalized image
BaseImageArray = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, BaseSeriesDescription, 'avgImageArray.mat')).avgImageArray;

Nz = size(ImageArray, 3);
BaseNz = size(BaseImageArray, 3);


%% First inspect images

figure
imshow(squeeze(BaseImageArray(120,:,:)), [])

figure
imshow(squeeze(ImageArray(5,:,:)), [])


%% Define slice range

slices = 10:12;


%% Define ROI

mask = zeros(size(ImageArray));

for sl = slices

    basesl = (BaseNz/Nz)*(2*sl-1)/2;

    f1=figure;
    tiledlayout(1,1);
    nexttile();    
    imshow(BaseImageArray(:,:,basesl), [])

    f2=figure;
    tiledlayout(1,1);
    nexttile();
    f2.Position = f2.Position + [580, 0, 0, 0];
    imshow(ImageArray(:,:, sl),[0 1])
    this_mask=drawfreehand();
    this_mask = createMask(this_mask);

    mask(:,:,sl) = this_mask;

    close(f1);
    close(f2);

end



%% Name and save ROI

% Name 
ROIname = 'UQ4M Lesion 3+4';

% Save
save(fullfile(projectfolder, 'Scripts', 'ROIs', [ROIname '.mat'] ), "mask");

