% MATLAB script to extra volume fraction information from ADC histograms

% Idea is to fit Gaussian distributions to histograms


%% Define image details and load ADC

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Output folder
OutputFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Outputs";

% Sample name
SampleName = '20241128_UQ1';

% Series description
SeriesDescription = '40u_DtiStandard_2012';

% Load ADC
ADC = load(fullfile(OutputFolder, SampleName, SeriesDescription, 'ADC.mat')).ADC;


%% Define ROI in each sample

SampleNumber = 1;
xs = 40:90;
ys = 34:84;
zs = 170:295;

% SampleNumber = 2;
% xs = 40:90;
% ys = 34:84;
% zs = ;


%% Histogram of ADC values

ADCvals = ADC(ys,xs,zs);

f = figure;
h = histogram(ADCvals(:), 400);
hold on
binedges = h.BinEdges;
bincentres = (1/2)*(binedges(1:end-1)+binedges(2:end));
counts = h.Values;



% Gaussian fitting to central peak
ADCpeak = mode(ADCvals, 'all');
ADCmin = 0.9e-3;
ADCmax = 1.2e-3;

bools = and(bincentres>ADCmin, bincentres<ADCmax);
thesecounts = counts(bools);
thesebincentres = bincentres(bools);

plot(thesebincentres, thesecounts, LineWidth=2);