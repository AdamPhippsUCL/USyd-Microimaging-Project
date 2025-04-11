% MATLAB script to test diffusion tensor calculation

clear;
projectfolder = pwd;

%% Image details

UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end

% Sample name
SampleName = '20250407_UQ5';

% Series description
SeriesDescription = '40u_DtiSE_2012_SPOIL10% (20 micron)';


%% Load image and dinfo

ImageArray = load(fullfile(ImagingDataFolder, SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;

%% Preprocessing

% b=0 slices
b0slices = [dinfo(:).DiffusionBValue] == 0;

% b0 image
b0img = ImageArray(:,:,:,b0slices);
b0img = mean(b0img, 4);

% b0 values (effective)
b0vals = [dinfo(b0slices).DiffusionEffBValue];
b0val = b0vals(1);

% b>0 slice 
bslices = [dinfo(:).DiffusionBValue] ~= 0;

% b>0 images
bimgs = ImageArray(:,:,:,bslices);

% b values (effective)
bvals = [dinfo(bslices).DiffusionEffBValue];
bvals = transpose(bvals);

% Relative effective b values
effbvals = bvals-b0val;

% Directions
direcs = [dinfo(bslices).DiffusionDirection];
direcs = reshape(direcs, [3,sum(bslices)]);
direcs = transpose(direcs);

% Normalise
norms = sqrt(sum(direcs.^2, 2));
direcs = direcs./norms;

% 
% % Select small image region for testing
% xs = 300:340;
% ys = 50:200;
% zs = 120;
% 
% b0img = b0img(zs,ys,xs);
% bimgs = bimgs(zs,ys,xs,:);
% 

%% Diffusion tensor and FA calculation

% Test DTI function
D_tensors = computeDiffusionTensor(bimgs, b0img, effbvals, direcs);

% Calculate FA
FA = computeFA(D_tensors);
FA = squeeze(FA);


%% Save results

Dx = real(squeeze(D_tensors(:,:,:,1)));
Dy = real(squeeze(D_tensors(:,:,:,2)));
Dz = real(squeeze(D_tensors(:,:,:,3)));
D = (1/3)*(Dx + Dy + Dz);

outfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', SeriesDescription);
mkdir(outfolder);

save(fullfile(outfolder, 'Dx.mat'), 'Dx');
save(fullfile(outfolder, 'Dy.mat'), 'Dy');
save(fullfile(outfolder, 'Dz.mat'), 'Dz');
save(fullfile(outfolder, 'D.mat'), 'D');
save(fullfile(outfolder, 'FA.mat'), 'FA');

% Diffusion weighted FA
dwFA = D.*FA;
save(fullfile(outfolder, 'dwFA.mat'), 'dwFA');
