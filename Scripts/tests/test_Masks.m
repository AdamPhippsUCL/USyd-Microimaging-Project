% Test script to look at T2/D/R measurements in epithelium and stroma

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

%% Parameter Map

% Sample name
samplename = '20241218_UQ3';

% Series description
seriesdescription = '3DMGE_20u (40 micron)';

% Model type
modeltype = 'T2';

% Parameter
parameter = 'T2';

% Load map
load(fullfile(projectfolder, 'Outputs', 'Model Fitting', samplename, modeltype, seriesdescription, [parameter '.mat'] ));

% Load voxel coordinates
highVoxelCoordinates = load(fullfile(projectfolder, 'Imaging Data', 'MAT', samplename, seriesdescription, 'axialVoxelCoordinates.mat')).VoxelCoordinates;

%% Mask

% Mask name
maskname = 'Test_Stroma1';

% Load mask
mask = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', samplename, modeltype, seriesdescription, [maskname '.mat'] )).img;
mask = double(mask);

% figure
% histogram(T2(logical(mask)), 50)
% xlim([0, 50])
% title(maskname)
% 


%% Low res parameter map

modeltype = 'RDI - 2 compartment - 4 param';
schemename = 'UQ3 Full';
fittingtechnique = 'MLP';

paramfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', samplename, modeltype, schemename, fittingtechnique);

% Load parameters
fIC = load(fullfile(paramfolder, 'fIC.mat')).fIC;
fEES = load(fullfile(paramfolder, 'fEES.mat')).fEES;
R = load(fullfile(paramfolder, 'R.mat')).R;
dIC = load(fullfile(paramfolder, 'dIC.mat')).dIC;
dEES = load(fullfile(paramfolder, 'dEES.mat')).dEES;




%% Upscale parameter maps

fIC_highres = imresize3(fIC, size(mask, 1:3));
save(fullfile(paramfolder, 'fIC_highres.mat'), 'fIC_highres')
fEES_highres = imresize3(fEES, size(mask, 1:3));
save(fullfile(paramfolder, 'fEES_highres.mat'), 'fEES_highres')
dIC_highres = imresize3(dIC, size(mask, 1:3));
save(fullfile(paramfolder, 'dIC_highres.mat'), 'dIC_highres')
dEES_highres = imresize3(dEES, size(mask, 1:3));
save(fullfile(paramfolder, 'dEES_highres.mat'), 'dEES_highres')
R_highres = imresize3(R, size(mask, 1:3));
save(fullfile(paramfolder, 'R_highres.mat'), 'R_highres')

meanfIC = mean(fIC_highres(logical(mask)))
meanR = mean(R_highres(logical(mask)))
meandIC = mean(dIC_highres(logical(mask)))



% % Doing it properly with voxel coordinates (there's a bug somewhere)

% lowVoxelCoordinates = load(fullfile(projectfolder, 'Imaging Data', 'MAT', samplename, lowseriesdescription, 'axialVoxelCoordinates.mat')).VoxelCoordinates;
% 
% % Adjusting to centre of voxels
% maxs = squeeze(lowVoxelCoordinates(end,end,end,:));
% mins = squeeze(lowVoxelCoordinates(1,1,1,:));
% 
% centrelowVoxelCoordinates = zeros(size(lowVoxelCoordinates));
% 
% centrelowVoxelCoordinates(:,:,:,1) = lowVoxelCoordinates(:,:,:,1) + (1/2)*(maxs(1)-mins(1))/(size(lowVoxelCoordinates,1)-1);
% centrelowVoxelCoordinates(:,:,:,2) = lowVoxelCoordinates(:,:,:,2) + (1/2)*(maxs(2)-mins(2))/(size(lowVoxelCoordinates,2)-1);
% centrelowVoxelCoordinates(:,:,:,3) = lowVoxelCoordinates(:,:,:,3) + (1/2)*(maxs(3)-mins(3))/(size(lowVoxelCoordinates,3)-1);

% % Resample mask at low coorindates
% lowmask = resample(mask, highVoxelCoordinates, lowVoxelCoordinates);
% lowmask = double(lowmask==1);

% % QUICK METHOD: reshaping mask
% lowmask = imresize3(mask, size(lowimg,1:3), method = "linear");
% 
% 
% 
% meanparam = mean(parammap(logical(lowmask)), 'all')


% %% Testing data quality
% 
% ImageFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN', samplename);
% 
% % Series descriptions
% seriesdescriptions = {
%     '40u_verdict_seq1_LR',...
%     '40u_verdict_seq2_LR',...
%     '40u_verdict_seq3_LR',...
%     '40u_verdict_seq4_LR',...
%     '40u_verdict_seq5_LR'};
% 
% figure
% 
% for sindx = 1:length(seriesdescriptions)
% 
%     sd = seriesdescriptions{sindx};
% 
%     % Load denoised image
%     img = load(fullfile(ImageFolder, sd, 'axialImageArray.mat')).ImageArray;
%     dinfo = load(fullfile(ImageFolder, sd, 'axialdinfo.mat')).dinfo;
% 
%     bvals = [dinfo(:).DiffusionBValue];
% 
%     % b0 imgs
%     b0bools = (bvals==0);
%     b0imgs = img(:,:,:,b0bools);
%     b0img = mean(b0imgs,4);
% 
%     b0signal = mean(b0img(logical(lowmask)));
% 
%     % b imgs
%     bbools = (bvals > 0);
%     bimgs = img(:,:,:,bbools);
%     bimg = mean(bimgs,4);
%     bval = bvals(find(bbools,1));
% 
%     bsignal = mean(bimg(logical(lowmask)));
% 
%     scatter(bval, bsignal/b0signal)
%     hold on
% 
% 
% end