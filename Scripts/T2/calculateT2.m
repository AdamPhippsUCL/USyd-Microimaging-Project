% MATLAB script to perform T2 calculation

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241213_WSU1';

% Series description
SeriesDescription = '160001';

%% Load image + info

ImageArray = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialImageArray.mat')).ImageArray;
dinfo = load(fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription, 'axialdinfo.mat')).dinfo;

% TE vector
TEvec = [dinfo(:).EchoTime];

%% T2 calculation

[T2, S0] = calcT2(ImageArray, TEvec);

% Clip unrealistic values
T2(T2<0)=0;
T2(T2>1000)=1000;
T2(isinf(T2))=0;
T2(isnan(T2))=0;



%% Save T2 image

% Output folder 
outputfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Outputs\Model Fitting";
folder = fullfile(outputfolder, SampleName,  'T2', SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'T2.mat'), 'T2'); 




%% Display example images

slice = 450;

figure;
imgslice = ImageArray(:,:,slice, 1);
imshow(imgslice, [])

figure;
T2slice = T2(:,:,slice);
imshow(T2slice, [0, 20]);


% %% Experiment with histograms
% 
% T2section = T2(80:160, 80:160, 400:500);
% 
% figure;
% binedges = linspace(0,20,501);
% h=histogram(T2section(:), BinEdges = binedges);
% hold on
% counts = h.Values;
% bincentres = (1/2)*(binedges(1:end-1)+binedges(2:end));
% binspacing = binedges(2)-binedges(1);
% 
% coeffs = fitGaussians( ...
%     counts, ...
%     bincentres, ...
%     N=1, ...
%     beta0guess = [0.8,8,2], ...
%     lb = [0,0,0], ...
%     ub = [1,20,20] );
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
% % test = stdfilt(T2section);
% % figure;
% % imshow(test, [0 2])
% % 
% % 
% % figure;
% % histogram(test(:), BinEdges = linspace(0,10,1001))
% % 
% % mask = (test>0.5);
% % 
% % figure;
% % imshow(double(mask), [])
