% MATLAB Script to create new training data for restricted diffusion imaging models

clear;

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";


% Number of training samples per parameter
N = 50;

% Noise 
noisetype = 'None';
sigma0 = 0.01;
T2 = 10000;

% Training data folder
savedata = true;
TrainingDataFolder = fullfile(projectfolder, 'Scripts', 'RDI', 'MLP', 'training data');


%% Protocol

% Model type
modeltype = 'RDI - 1 compartment - 2 param';

% Scheme name
schemename = '20250224_UQ4 AllDELTA';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";

%% Parameter ranges


switch modeltype

    case 'RDI - 1 compartment - 2 param'

        Nparam = 2;

        Rmin = 1;
        Rmax = 50;

        dICmin = 0.1;
        dICmax = 3;

        minparams = [Rmin, dICmin];
        maxparams = [Rmax, dICmax];


    case 'RDI - 2 compartment - 3 param'

        Nparam = 3;

        fICmin = 0;
        fICmax = 1;

        Rmin = 1;
        Rmax = 50;

        dmin = 0.1;
        dmax = 3;

        minparams = [fICmin, Rmin, dmin];
        maxparams = [fICmax, Rmax, dmax];


    case 'RDI - 2 compartment - 4 param'

        Nparam = 4;

        fICmin = 0;
        fICmax = 1;

        Rmin = 5;
        Rmax = 50;

        dICmin = 0.2;
        dICmax = 3;

        dEESmin = 0.2;
        dEESmax = 3;

        minparams = [fICmin, Rmin, dICmin, dEESmin];
        maxparams = [fICmax, Rmax, dICmax, dEESmax];

end

% AMICO
Rvals = linspace(0.1,20.1,21);


%% Create training data

Nvoxel = N^Nparam;

outputfolder = fullfile(TrainingDataFolder, modeltype, schemename);

createRDIdata( ...
    modeltype,...
    schemename,...
    Nvoxel = Nvoxel,...
    noisetype = noisetype,...
    sigma0 = sigma0,...
    T2 = T2,...
    minparams = minparams,...
    maxparams = maxparams,...
    Rvals = Rvals,... 
    schemesfolder=schemesfolder,...
    savedata=savedata,...
    outputfolder=outputfolder...
    );
