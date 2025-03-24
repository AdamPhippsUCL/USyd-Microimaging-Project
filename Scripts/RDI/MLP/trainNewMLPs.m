% Script to train MLP models

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

TrainingDataFolder = fullfile(projectfolder, 'Scripts', 'RDI', 'MLP', 'training data');
ModelFolder = fullfile(projectfolder, 'Scripts', 'RDI', 'MLP', 'models');

% Specify protocols 
modeltypes = {'RDI - 1 compartment - 2 param'};
schemenames = {'20250224_UQ4 ShortDELTA'};

% scale params
scaleparams = true;

% Run MLP training
for indx = 1:length(modeltypes)
    
    modeltype = modeltypes{indx};
    schemename = schemenames{indx};

    datafolder = fullfile(TrainingDataFolder, modeltype, schemename);
    modelfolder = fullfile(ModelFolder, modeltype, schemename);

    trainMLP( ...
        datafolder,...
        modelfolder,...
        maxepochs = 150, ...
        scaleparams = scaleparams...
        );
    
end