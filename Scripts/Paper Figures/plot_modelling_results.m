% Script to compare results from direct modelling, to predicted model
% parameters from ESL component fractions

clear;
projectfolder = pwd;

%%

SampleNames = {...
    '20250224_UQ4',...
    '20250407_UQ5',...
    '20250414_UQ6',...
    '20250522_UQ7',...
    '20250523_UQ8',...
    '20250524_UQ9'...
};

ModelName = 'RDI - 2 compartment - 4 param (S0)';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

% Load ESL modelling estimates
ESL_Model_RESULTS = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'Modelling', 'RESULTS.mat')).RESULTS;
S_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'S'), strcmp({ESL_Model_RESULTS(:).ModelType}, ModelName))).ModelParams;
E_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'G'), strcmp({ESL_Model_RESULTS(:).ModelType}, ModelName))).ModelParams;
L_params = [1, 2]; % ADC parameters...




% Load model parameters from direct fitting, and composition array
for sampleindx = 1:length(SampleNames)

    SampleName = SampleNames{sampleindx};

    % Load composition array
    sample_composition = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;
  
    switch ModelName
        
        case 'RDI - 2 compartment - 4 param (S0)'

            outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, ModelName, schemename, fittingtechnique);
            fit_fIC = load(fullfile(outputfolder, 'fIC.mat')).fIC;
            fit_dIC = load(fullfile(outputfolder, 'dIC.mat')).dIC;
            fit_R = load(fullfile(outputfolder, 'R.mat')).R;
            fit_dEES = load(fullfile(outputfolder, 'dEES.mat')).dEES;

            % Predicted model parameters
            pred_fIC = S_params(1).*sample_composition(:,:,:,1) + E_params(1).*sample_composition(:,:,:,2);
            pred_R = S_params(2).*sample_composition(:,:,:,1) + E_params(2).*sample_composition(:,:,:,2);
            pred_dIC = S_params(3).*sample_composition(:,:,:,1) + E_params(3).*sample_composition(:,:,:,2);
            pred_dEES = S_params(4).*sample_composition(:,:,:,1) + E_params(4).*sample_composition(:,:,:,2) + L_params(2).*sample_composition(:,:,:,3);

    end

    % Inspect differences in pred and fit maps here...

    



    this_composition = reshape(sample_composition,[prod(size(sample_composition, 1:3)), 3]);
    bool =  sum(this_composition,2)>0;
    this_composition = this_composition(bool, :);
    bool = reshape(bool, size(sample_composition, 1:3));


    switch ModelName

        case 'RDI - 2 compartment - 4 param (S0)'

        % Extract modelling and prediction results for non-zero composition
        this_fit_fIC = fit_fIC(bool);
        this_fit_R = fit_R(bool);
        this_fit_dIC = fit_dIC(bool);
        this_fit_dEES = fit_dEES(bool);
    
        this_pred_fIC = pred_fIC(bool);
        this_pred_R = pred_R(bool);
        this_pred_dIC = pred_dIC(bool);
        this_pred_dEES = pred_dEES(bool);
    
        % Append to parameter array
        if sampleindx == 1       
            composition = this_composition;
            pred_params = [this_pred_fIC'; this_pred_R'; this_pred_dIC'; this_pred_dEES'];
            fit_params = [this_fit_fIC'; this_fit_R'; this_fit_dIC'; this_fit_dEES'];
        else
            composition = cat(1, composition, this_composition);
            pred_params = cat(2, pred_params, [this_pred_fIC'; this_pred_R'; this_pred_dIC'; this_pred_dEES']);
            fit_params = cat(2, fit_params, [this_fit_fIC'; this_fit_R'; this_fit_dIC'; this_fit_dEES']);
        end
        
    end


end



%% Display results

% Reorganise composition array
new_composition = zeros(size(composition));
new_composition(:,1) = composition(:,2);
new_composition(:,2) = composition(:,1);
new_composition(:,3) = composition(:,3);


figure
scatter(fit_params(1,:), pred_params(1,:), CData=new_composition)

figure
scatter(fit_params(2,:), pred_params(2,:), CData=new_composition)

figure
scatter(fit_params(3,:), pred_params(3,:), CData=new_composition)

figure
scatter(fit_params(4,:), pred_params(4,:), CData=new_composition)