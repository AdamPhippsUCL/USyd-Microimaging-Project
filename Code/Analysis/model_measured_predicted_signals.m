% Script to apply diffusion signal model to measured and predicted signals

clear;
projectfolder = pwd;

%% Sample and image details

% % Sample groups
% Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };
% Cancer_3 = {'4B', '4M'};
% Cancer_4 = {'6N'};
% 
% % Set group interested in
% group = 'Cancer_3';

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

% switch group
%     case 'Benign'
%         Bools = ismember(SampleNums, Benign);
%     case 'Cancer_3'
%         Bools = ismember(SampleNums, Cancer_3);
%     case 'Cancer_4'
%         Bools = ismember(SampleNums, Cancer_4);
% end
% Nvoxel = sum(Bools==1);


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

    % % Get signals from group
    % this_measured = this_measured(Bools);
    % this_pred = this_pred(Bools);

    MeasuredSignals(:,seriesindx) = this_measured;
    PredictedSignals(:,seriesindx) = this_pred;

end

% COMP = COMP(Bools, :);


%% Run modelling

% BALL+SPHERE Model

modelname = 'Ball+Sphere';
fittingtechnique = 'LSQ';
Nparam = 5;

fs = 0.5; fslb = 0; fsub = 1;
R = 6.5; Rlb = 6.4; Rub = 6.6;
Ds = 0.55; Dslb = 0.54; Dsub = 0.56;
Db = 1; Dblb = 0.1; Dbub = 3;
S0 = 1; S0lb = 0.9; S0ub = 1.1;

beta0 = [fs,R,Ds,Db,S0];
lb = [fslb,Rlb,Dslb,Dblb,S0lb];
ub = [fsub,Rub,Dsub,Dbub,S0ub];

% Regularisation
lambda0=0e-2;
lambda = lambda0*ones(1,Nparam);


% == MEASURED

Y = reshape(MeasuredSignals, [Nvoxel, 1 , nscheme]);

[measured_fs, ~, ~, measured_Db, ~, ~] = diffusion_model_fit( ...
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


% == PREDICTED

Y = reshape(PredictedSignals, [Nvoxel, 1 , nscheme]);

[pred_fs, ~, ~, pred_Db, ~, ~] = diffusion_model_fit( ...
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