% Script to train MLP models

% Base folder
basef = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Scripts\VERDICT");

% Specify protocols 
modeltypes = {'Original VERDICT'};
schemenames = {'UQ Scheme v2'};

TrainingDataFolder = [basef '/MLP/training data'];
ModelFolder = [basef '/MLP/models'];

for indx = 1:length(modeltypes)
    
    modeltype = modeltypes{indx};
    schemename = schemenames{indx};

    datafolder = [TrainingDataFolder '/' modeltype '/' schemename ];
    modelfolder = [ModelFolder '/' modeltype '/' schemename ];

    trainMLP( ...
        datafolder,...
        modelfolder...
        );
    
end