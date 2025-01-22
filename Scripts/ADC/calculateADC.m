% MATLAB Script to perform ADC calculation

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241218_UQ3';

% Series description
SeriesDescription = '40u_DtiStandard_2012';


%% Load image + info

ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;


%% Reformat images 

% Average b=0 images
whereb0 = ([dinfo(:).DiffusionBValue] == 0);
b0imgs = ImageArray(:,:,:,whereb0);
b0img = mean(b0imgs, 4);

% Average b>0 images
whereb = ([dinfo(:).DiffusionBValue] ~= 0);
bimgs = ImageArray(:,:,:,whereb);
bimg = mean(bimgs, 4);
bval = dinfo(find(whereb,1)).DiffusionBValue;

% Gaussian smoothing
gaussfilt = false;
sigma = 0.5;
if gaussfilt
    b0img = imgaussfilt(b0img);
    bimg = imgaussfilt(bimg);
end

% Input image
inputimg = cat(4, b0img, bimg); 
bvec = [40 bval];


%% ADC calculation

[ADC, S0] = calcADC(inputimg,  bvec);


%% Save ADC image

% Output folder 
outputfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Outputs\Model Fitting";
folder = fullfile(outputfolder, SampleName,  'ADC', SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'ADC.mat'), 'ADC'); 

% Meta
Meta = struct();
Meta.gaussfilt = gaussfilt;
Meta.sigma = sigma;
save(fullfile(folder, 'ADC Meta.mat'), 'Meta'); 


%% Example images
slice = 320;
figure;
b0slice = b0img(:,:,slice);
imshow(b0slice, [0 prctile(b0slice(:), 99)]);
figure;
bslice = bimg(:,:,slice);
imshow(bslice, [0 prctile(bslice(:), 99)]);
figure;
ADCslice = ADC(:,:,slice);
imshow(ADCslice, [0 prctile(ADCslice(:), 99)]);
title('ADC')
% colorbar;


% %% Experimenting with histograms
% 
% ADCvals = ADC(40:85, 40:85,200:300);
% figure;
% h=histogram(ADCvals(:), 250);
% hold on
% binedges = h.BinEdges;
% counts = h.Values;
% bincentres = (1/2)*(binedges(1:end-1)+binedges(2:end));
% binspacing = binedges(2)-binedges(1);
% 
% 
% coeffs = fitGaussians( ...
%     counts, ...
%     bincentres, ...
%     N=2, ...
%     beta0guess = [0.5,1e-3,5e-4, 0.5, 1.6e-3,5e-4], ...
%     lb = [0,0,0,0,0,0,0,0,0], ...
%     ub = [1,20,20,1,20,20, 1,20,20] );
% 
% 
% % Normalize weights
% coeffs(1:3:end) = coeffs(1:3:end)/sum(coeffs(1:3:end));
% 
% % Create distribution pdfs
% N = length(coeffs)/3;
% dists = zeros(N, length(bincentres));
% for n = 1:N
%     dists(n,:) = coeffs(3*n-2)*normpdf(bincentres, coeffs(3*n-1), coeffs(3*n));
%     plot(bincentres, binspacing*sum(counts(:))*dists(n,:));
% end
% 
% plot(bincentres, binspacing*sum(dists,1)*sum(counts(:)), LineWidth = 2)
% 
% % mask = (ADCslice<1.1e-3);
% % 
% % figure
% % imshow(mask, [])
