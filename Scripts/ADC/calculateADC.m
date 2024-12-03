% MATLAB Script to perform ADC calculation

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241128_UQ1';

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
bvec = [0 bval];


%% ADC calculation

[ADC, S0] = calcADC(inputimg,  bvec);


%% Save ADC image

% Output folder 
outputfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Outputs";
folder = fullfile(outputfolder, SampleName, SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'ADC.mat'), 'ADC'); 

% Meta
Meta = struct();
Meta.gaussfilt = gaussfilt;
Meta.sigma = sigma;
save(fullfile(folder, 'ADC Meta.mat'), 'Meta'); 


%% Example images
slice = 222;
figure;
b0slice = b0img(:,:,slice);
imshow(b0slice, [0 prctile(b0slice(:), 99)]);
figure;
bslice = bimg(:,:,slice);
imshow(bslice, [0 prctile(bslice(:), 99)]);
figure;
ADCslice = ADC(:,:,slice);
imshow(ADCslice, [0 prctile(ADCslice(:), 99)]);
% 
% 
% %% Experimenting with histograms
% 
% ADCvals = ADC(40:85, 40:85,200:300);
% figure;
% h=histogram(ADCvals(:), 250);
% binedges = h.BinEdges;
% bincentres = (1/2)*(binedges(1:end-1) + binedges(2:end));
% counts = h.Values;
% % 
% mask = (ADCslice<1.1e-3);
% 
% figure
% imshow(mask, [])
