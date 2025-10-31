% Script to apply diffusion signal model to measured and predicted signals

clear;
projectfolder = pwd;

%% Image details

% Images
SeriesDescriptions = {
    'SE_b0_SPOIL5% (DS)',...
    'STEAM_ShortDELTA_15 (DS)',...
    'STEAM_ShortDELTA_20 (DS)',...
    'STEAM_ShortDELTA_30 (DS)',...
    'STEAM_ShortDELTA_40 (DS)',...
    'STEAM_ShortDELTA_50 (DS)',...
    'STEAM_LongDELTA_40 (DS)',...
    'STEAM_LongDELTA_60 (DS)',...
    'STEAM_LongDELTA_80 (DS)',...
    'STEAM_LongDELTA_100 (DS)',...
    'STEAM_LongDELTA_120 (DS)'...
};
scheme = load(fullfile(projectfolder, "Schemes", "20250224_UQ4 AllDELTA.mat")).scheme;
nscheme = length(scheme);

folder =  fullfile(projectfolder, 'Outputs', 'Signals');
COMP = load(fullfile(folder, "COMP.mat")).COMP;
SampleNums = load(fullfile(folder, "SampleNums.mat")).SampleNums;
Nvoxel = length(SampleNums);

%% Load predicted and measured signals

% Initialise array for predicted signals
MeasuredSignals = ones(Nvoxel, nscheme);
PredictedSignals = ones(Nvoxel, nscheme);

for seriesindx = 2:length(SeriesDescriptions)
    
    SeriesDescription = SeriesDescriptions{seriesindx};

    bval = scheme(seriesindx).bval;
    DELTA = scheme(seriesindx).DELTA;

    % Load measured and predicted signals
    this_measured = load(fullfile(folder, SeriesDescription, "Measured.mat")).Measured;
    this_pred = load(fullfile(folder, SeriesDescription, "Predicted.mat")).Predicted;

    MeasuredSignals(:,seriesindx) = this_measured;
    PredictedSignals(:,seriesindx) = this_pred;

end


%% Run ADC modelling


% ============= ADC model

modelname = 'ADC';
fittingtechnique = 'LSQ';
Nparam = 2;

D = 1; Dlb = 0.1; Dub = 3;
S0 = 1; S0lb = 0.9; S0ub = 1.1;

beta0 = [S0, D];
lb = [S0lb, Dlb];
ub = [S0ub, Dub];

% Regularisation
lambda0=0e-2;
lambda = lambda0*ones(1,Nparam);


% == MEASURED

Y = reshape(MeasuredSignals, [Nvoxel, 1 , nscheme]);

[~, measured_D] = diffusion_model_fit( ...
    Y, ...
    scheme, ...
    modelname = modelname,...
    fittingtechnique = fittingtechnique,...
    Nparam = Nparam,...
    beta0=beta0,...
    lambda=lambda,...
    lb=lb,...
    ub=ub...
    );

% Save
folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Measured', modelname);
mkdir(folder);
save(fullfile(folder, 'D.mat'), 'measured_D')

% == PREDICTED

Y = reshape(PredictedSignals, [Nvoxel, 1 , nscheme]);

% Normalise for S0 
Y(:,:,2:end) = Y(:,:,2:end)./(sum(COMP,2));

[~, pred_D] = diffusion_model_fit( ...
    Y, ...
    scheme, ...
    modelname = modelname,...
    fittingtechnique = fittingtechnique,...
    Nparam = Nparam,...
    beta0=beta0,...
    lambda=lambda,...
    lb=lb,...
    ub=ub...
    );

% Save
folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Predicted', modelname);
mkdir(folder);
save(fullfile(folder, 'D.mat'), 'pred_D')




%% Run Ball+Sphere Modelling

modelname = 'Ball+Sphere';
fittingtechnique = 'LSQ';

Nparam = 5;
beta0 = [0.24, 6.4, 0.6, 0.8, 1];
lb = [0, 2, 0.4, 0.2, 1];
ub = [1, 12, 0.8, 3, 1];


% Regularisation (Don't regluarlise parameters which would change
% significantly in purely fluid voxel i.e. Db and fs)
lambda0 = 2e-2;
lambda = lambda0*[0,1,1,0,1];


% == PREDICTED

Y = reshape(PredictedSignals, [Nvoxel, 1 , nscheme]);

% Normalise for S0
Y(:,:,2:end) = Y(:,:,2:end)./(sum(COMP,2));

[pred_fs, pred_R, pred_Ds, pred_Db, ~, ~] = diffusion_model_fit( ...
    Y, ...
    scheme, ...
    modelname = modelname,...
    fittingtechnique = fittingtechnique,...
    Nparam = Nparam,...
    beta0=beta0,...
    lambda=lambda,...
    lb=lb,...
    ub=ub...
    );

% Save
folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Predicted', modelname);
mkdir(folder);
save(fullfile(folder, 'fs.mat'), 'pred_fs')
save(fullfile(folder, 'Db.mat'), 'pred_Db')
save(fullfile(folder, 'Ds.mat'), 'pred_Ds')
save(fullfile(folder, 'R.mat'), 'pred_R')


% == MEASURED

Y = reshape(MeasuredSignals, [Nvoxel, 1 , nscheme]);

[measured_fs, measured_R, measured_Ds, measured_Db, ~, ~] = diffusion_model_fit( ...
    Y, ...
    scheme, ...
    modelname = modelname,...
    fittingtechnique = fittingtechnique,...
    Nparam = Nparam,...
    beta0=beta0,...
    lambda=lambda,...
    lb=lb,...
    ub=ub...
    );

% Save
folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'Measured', modelname);
mkdir(folder);
save(fullfile(folder, 'fs.mat'), 'measured_fs')
save(fullfile(folder, 'Db.mat'), 'measured_Db')
save(fullfile(folder, 'Ds.mat'), 'measured_Ds')
save(fullfile(folder, 'R.mat'), 'measured_R')


