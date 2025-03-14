% MATLAB Script to create new training data for restricted diffusion imaging models

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
TrainingDataFolder = fullfile(projectfolder, 'Scripts', 'RDI', 'MLP', 'training data');


%% Protocol

% Model type
modeltype = 'RDI - 2 compartment - 4 param';

% Scheme name
schemename = '20250224_UQ4 MixedDELTA2';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";

%% Parameter ranges

% Intracellular compartment
fICs = [0, 1];
dICs = [0.5, 3];
Rs = [4, 16];
% Extracellular compartment
dEESs = [0.5, 3];

% AMICO
Rvals = linspace(0.1,20.1,21);
muRs = [4, 12];
sigmaRs = [1, 2];




%% Create training data

outputfolder = fullfile(TrainingDataFolder, modeltype, schemename);

createRDIdata( ...
    modeltype,...
    schemename,...
    Nvoxel = Nvoxel,...
    noisetype = noisetype,...
    sigma0 = sigma0,...
    T2 = T2,...
    fICs = fICs,...
    dICs = dICs,...
    Rs = Rs,...
    Rvals = Rvals,...
    muRs = muRs,...
    sigmaRs = sigmaRs,... 
    dEESs = dEESs,...
    schemesfolder=schemesfolder,...
    savedata=savedata,...
    outputfolder=outputfolder...
    );
