% Script to create new training data for MLP model training

% Base folder
basef = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Scripts\VERDICT");

% Number of training data samples
Nvoxel = 20000;

% Noise 
noisetype = 'Rice';
sigma0 = 0.01;
T2 = 10000;

% Training data folder
savedata = true;
TrainingDataFolder = [basef '/MLP/training data'];

%% Protocols

% If multiple, define as a cell array of char arrays
modeltypes = {'Original VERDICT'};
schemenames = {'UQ Scheme v2'};

schemesfolder = [basef '/Schemes'];



%% Cell sizes

% Cell Radius distribution R~Normal( muR, sigmaR)
muRmin = 6;
muRmax = 9;
sigmaRmin = 1.5;
sigmaRmax = 2.5;


%% Create training data

for indx = 1:length(modeltypes)

    modeltype = modeltypes{indx};
    schemename = schemenames{indx};

    disp([modeltype ' ' schemename])

    outputfolder = [TrainingDataFolder '/' modeltype '/' schemename '/'];
    
    createVERDICTdata( ...
        modeltype,...
        schemename,...
        Nvoxel = Nvoxel,...
        noisetype = noisetype,...
        sigma0 = sigma0,...
        T2 = T2,...
        randmuRs=[muRmin, muRmax],...
        randsigmaRs=[sigmaRmin, sigmaRmax],...
        schemesfolder=schemesfolder,...
        savedata=savedata,...
        outputfolder=outputfolder...
        );

end

