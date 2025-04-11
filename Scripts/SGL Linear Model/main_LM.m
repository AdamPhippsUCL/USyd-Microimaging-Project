% Script to apply linear model to model parameter maps

clear;
projectfolder = pwd;


%% Parameter map details

samplename = '20250224_UQ4';

% modeltype = 'ADC';
modeltype = 'DKI';

schemename = '20250224_UQ4 AllDELTA';
scheme = load(fullfile(projectfolder, 'Schemes', [schemename '.mat'])).scheme;
nscheme = length(scheme);

fittingtechnique = 'LSQ';

parameter = 'K';

% Parameter map
paramfolder = fullfile(projectfolder, "Outputs", "Model Fitting", samplename, modeltype, schemename, fittingtechnique);
parammap = load(fullfile(paramfolder,[parameter '.mat'])).(parameter);
szmap = size(parammap);
parammap_flat = reshape(parammap, [prod(szmap), 1]);


% Intercept in model
switch parameter
    case {'fIC', 'C'}
        intercept = false;
    case {'ADC','D', 'R', 'dEES', 'dIC', 'd', 'K', 'FA'}
        intercept = true;
end



%% Load masks

seriesdescription = '3DMGE_20u';
maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, seriesdescription);
GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;
szmask = size(GLANDULAR);


%% Sample mask

switch samplename

    case '20250224_UQ4'

        % Cylinder centred at (128, 114)      
        samplemask = zeros(szmask);
        
        [Xs, Ys] = meshgrid( ...
           1:szmask(2), ...
           1:szmask(1) ...
            );
        
        samplemask(:,:,:) = repmat( ( (Xs-128).^2 + (Ys-114).^2 < 75^2 ), 1, 1, szmask(3));

end


%% Construct composition image

COMPOSITION = zeros([szmap, 3]);
ResFactor = szmask./szmap;
Nvoxel = prod(ResFactor);

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
        for slindx = 1:szmap(3)

            rows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
            cols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
            slices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));

            % Test in sample
            if ~all(logical(samplemask(rows, cols, slices)), "all")
                continue
            end

            % COMPOSITION         
            COMPOSITION(rindx, cindx, slindx, 1) = sum(double(STROMA(rows, cols, slices)), "all")/Nvoxel;
            COMPOSITION(rindx, cindx, slindx, 2) = sum(double(GLANDULAR(rows, cols, slices)), "all")/Nvoxel;
            COMPOSITION(rindx, cindx, slindx, 3) = sum(double(LUMEN(rows, cols, slices)), "all")/Nvoxel;        

        end
    end
end

% Select voxels with non-zero composition
composition = reshape(COMPOSITION,[prod(szmap), 3]);
bool =  sum(composition,2)>0;
composition = composition(bool, :);


%% Linear model

% Predictors and labels
X = composition;
paramvals = parammap_flat(bool);
y = paramvals;

% Fit a multiple linear regression model
mdl = fitlm(X, y,'Intercept', intercept);

y_pred = predict(mdl, X); % Get predicted values
r = corr(y, y_pred);

figure;
scatter(y, y_pred, 'filled'); % Plot actual vs predicted
xlabel('Actual Values');
ylabel('Predicted Values');
title('Actual vs. Predicted Values (LM)');
hold on;
plot(y, y, '--'); % 45-degree reference line


% Standardised coefficients
switch intercept
    case true
        b = mdl.Coefficients.Estimate(2:end); % Exclude intercept

        sx = std(X);  % Standard deviation of predictors
        sy = std(y);  % Standard deviation of response
        beta_std = (b .* sx') / sy; % Standardised coefficients
        disp(table(mdl.CoefficientNames(2:end)', beta_std, 'VariableNames', {'Predictor', 'Standardised_Coeff'}));
        
        % Standardise confidence intervals
        ci_raw = coefCI(mdl, 0.05);  % 95% CI (default)
        ci_std = (ci_raw(2:end, :) .* sx') / sy;  % Exclude intercept row
        % Display results in a table
        Standardised_Coeff_Table = table(mdl.CoefficientNames(2:end)', beta_std, ci_std(:,1), ci_std(:,2), ...
            'VariableNames', {'Predictor', 'Standardised_Coeff', 'Lower_CI', 'Upper_CI'});
        disp(Standardised_Coeff_Table)

    case false
        b = mdl.Coefficients.Estimate(1:end); % 
        
        sx = std(X);  % Standard deviation of predictors
        sy = std(y);  % Standard deviation of response
        beta_std = (b .* sx') / sy; % Standardised coefficients
        disp(table(mdl.CoefficientNames(1:end)', beta_std, 'VariableNames', {'Predictor', 'Standardised_Coeff'}));

        % Standardise confidence intervals
        ci_raw = coefCI(mdl, 0.05);  % 95% CI (default)
        ci_std = (ci_raw(1:end, :) .* sx') / sy;  
        Standardised_Coeff_Table = table(mdl.CoefficientNames(1:end)', beta_std, ci_std(:,1), ci_std(:,2), ...
            'VariableNames', {'Predictor', 'Standardised_Coeff', 'Lower_CI', 'Upper_CI'});
        disp(Standardised_Coeff_Table)       
end


% Results structure
RESULTS = struct();
RESULTS.modeltype = modeltype;
RESULTS.parameter = parameter;
RESULTS.mdl = mdl;

figure
bar(beta_std)
hold on
xticklabels({'S', 'G', 'L'})
% ylim([-1,1])
title([modeltype '; ' parameter])
errorbar(1:3, mean(ci_std,2), range(ci_std,2)/2, 'k', 'LineStyle', 'none', 'Marker', 'none')
ylabel('Standardised coefficient')



%% AIC

% Load resnorm
resnorm = load(fullfile(paramfolder, 'RESNORM.mat')).RESNORM;

% Number of parameters
switch modeltype
    case {'RDI - 2 compartment - 4 param'}
        Nparam = 4;
    case {'RDI - 2 compartment - 3 param', 'DKI'}
        Nparam = 3;
    case {'RDI - 1 compartment - 2 param', 'ADC'}
        Nparam = 2;
end

% AIC
AIC = (nscheme)*log((resnorm)/nscheme)-2*Nparam;
AICflat = reshape(AIC, [prod(szmap), 1]);
AICvals = AICflat(bool);

figure
histogram(AICvals)
