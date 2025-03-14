% Script to create new training data for MLP model training

% Base folder
basef = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Scripts\VERDICT");

% Number of training data samples
Nvoxel = 50000;

% Noise 
noisetype = 'Rice';
sigma0 = 0.02;
T2 = 10000;

% Training data folder
savedata = true;
TrainingDataFolder = [basef '/MLP/training data'];

%% Protocols

% If multiple, define as a cell array of char arrays
modeltypes = {'No VASC VERDICT (AMICO)'};
schemenames = {'20250224_UQ4 ShortDELTA'};

schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";


%% VERDICT model parameter ranges

fICs = [0, 1];
fVASCs = [0, 0];


% Cell Radius distribution R~Normal( muR, sigmaR)

Rvals = linspace(0.1, 15.1, 16);

muRmin = 3;
muRmax = 15;
Rs = [muRmin, muRmax];

sigmaRmin = 1;
sigmaRmax = 2;
sigmaRs = [sigmaRmin, sigmaRmax];

% Diffusivities
dIC=1;
dEES=1;
dVASC=8;


%% Create training data

for indx = 1:length(modeltypes)

    modeltype = modeltypes{indx};
    schemename = schemenames{indx};

    disp([modeltype ' ' schemename])

    outputfolder = fullfile(TrainingDataFolder, modeltype, schemename);
    
    createVERDICTdata( ...
        modeltype,...
        schemename,...
        Nvoxel = Nvoxel,...
        noisetype = noisetype,...
        sigma0 = sigma0,...
        T2 = T2,...
        fICs = fICs,...
        fVASCs = fVASCs,...
        Rs = Rs, ...
        Rvals = Rvals,....
        sigmaRs=sigmaRs,...
        dIC=dIC,...
        dEES=dEES,...
        dVASC=dVASC,...
        schemesfolder=schemesfolder,...
        savedata=savedata,...
        outputfolder=outputfolder...
        );

end

