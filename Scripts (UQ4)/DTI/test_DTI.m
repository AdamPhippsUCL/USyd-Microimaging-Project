% MATLAB script to test diffusion tensor calculation

clear;

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series description
SeriesDescription = '40u_DtiSE_2012_SPOIL10% (20 micron)';

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


%% Preprocessing

% b=0 slices
b0slices = [dinfo(:).DiffusionBValue] == 0;

% b0 image
b0img = ImageArray(:,:,:,b0slices);
b0img = mean(b0img, 4);

% b>0 slice 
bslices = [dinfo(:).DiffusionBValue] ~= 0;

% b>0 images
bimgs = ImageArray(:,:,:,bslices);

% b values
bvals = [dinfo(bslices).DiffusionEffBValue];
bvals = transpose(bvals);

% Directions
direcs = [dinfo(bslices).DiffusionDirection];
direcs = reshape(direcs, [3,sum(bslices)]);
direcs = transpose(direcs);
% Normalise
norms = sqrt(sum(direcs.^2, 2));
direcs = direcs./norms;


% % SELECT small image region
% xs = 1:640;
% ys = 50:200;
% zs = 120;
% 
% b0img = b0img(zs,ys,xs);
% bimgs = bimgs(zs,ys,xs,:);


%% Diffusion tensor and FA calculation

% Test DTI function
D_tensors = computeDiffusionTensor(bimgs, b0img, bvals, direcs);

% Calculate FA
FA = computeFA(D_tensors);
FA = squeeze(FA);


%% Save Dx, Dy, Dz, and FA

Dx = real(squeeze(D_tensors(:,:,:,1)));
Dy = real(squeeze(D_tensors(:,:,:,2)));
Dz = real(squeeze(D_tensors(:,:,:,3)));
D = (1/3)*(Dx + Dy + Dz);

projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

outfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', SeriesDescription);
mkdir(outfolder);

save(fullfile(outfolder, 'Dx.mat'), 'Dx');
save(fullfile(outfolder, 'Dy.mat'), 'Dy');
save(fullfile(outfolder, 'Dz.mat'), 'Dz');
save(fullfile(outfolder, 'D.mat'), 'D');
save(fullfile(outfolder, 'FA.mat'), 'FA');


%% Diffusion weighted FA

dwFA = D.*FA;

save(fullfile(outfolder, 'dwFA.mat'), 'dwFA');