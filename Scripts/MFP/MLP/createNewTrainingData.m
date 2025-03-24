% Script to create training data for mean free path imaging

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

% MFP folder
mfpfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\MFP-Processing");

% Number of training samples
Nvoxel = 50000;

% Noise 
noisetype = 'Rice';
sigma0 = 0.01;
T2 = 10000;

% Training data folder
savedata = true;
TrainingDataFolder = fullfile(projectfolder, 'Scripts', 'MFP',  'MLP', 'training data');


%% Protocols

% If multiple, define as a cell array of char arrays
modeltypes = {'MFP v2'};
schemenames = {'20250224_UQ4 ShortDELTA'};

% Set schemes folder
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";

%% Training data settings

% Range of MFP
Rmin = 4;
Rmax = 40;

% Diffusivity range
Dmin = 0.5;
Dmax = 3;


%% Create training data

for indx = 1:length(modeltypes)

    modeltype = modeltypes{indx};
    schemename = schemenames{indx};

    disp([modeltype ' ' schemename])

    outputfolder = fullfile(TrainingDataFolder, modeltype, schemename);
    
    createMFPdata( ...
        modeltype,...
        schemename,...
        Nvoxel = Nvoxel,...
        noisetype = noisetype,...
        sigma0 = sigma0,...
        T2 = T2,...
        Rs = [Rmin, Rmax],...
        Ds = [Dmin, Dmax],...
        schemesfolder=schemesfolder,...
        savedata=savedata,...
        outputfolder=outputfolder...
        );

end
