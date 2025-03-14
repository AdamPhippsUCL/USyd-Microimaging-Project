% MATLAB Script to create new training data for ADC model

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";


% Number of training samples
Nvoxel = 50000;

% Noise 
noisetype = 'Rice';
sigma0 = 0.02;
T2 = 10000;

% Training data folder
savedata = true;
TrainingDataFolder = fullfile(projectfolder, 'Scripts', 'ADC', 'MLP', 'training data');


%% Protocol

% Model type
modeltype = 'ADC';

% Scheme name
schemename = '20250224_UQ4 ShortDELTA';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";

%% Parameter ranges

% S0
S0 = [0.9, 1.1];
ADC = [0 4e-3];

