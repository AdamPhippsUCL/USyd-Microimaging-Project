% Script to train MLP models

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

% MFP folder
mfpfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\MFP-Processing");

% Specify protocols 
modeltypes = {'MFP v1'};
schemenames = {'UQ3 Full'};

TrainingDataFolder = fullfile(projectfolder, 'Scripts', 'MFP',  'MLP', 'training data');
ModelFolder = fullfile(projectfolder, 'Scripts', 'MFP',  'MLP', 'models');

for indx = 1:length(modeltypes)
    
    modeltype = modeltypes{indx};
    schemename = schemenames{indx};

    datafolder = fullfile(TrainingDataFolder, modeltype, schemename);
    modelfolder = fullfile(ModelFolder, modeltype, schemename);

    trainMLP( ...
        datafolder,...
        modelfolder,...
        maxepochs = 500 ...
        );
    
end