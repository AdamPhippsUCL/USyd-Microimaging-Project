% Script to apply linear model to model parameter maps

clear;
projectfolder = pwd;


%% Parameter map details

SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'}; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

modeltype = 'RDI - 2 compartment - 4 param (S0)'; % 'ADC', 'DKI', 'RDI - 2 compartment - 4 param (S0)'

schemename = '20250224_UQ4 AllDELTA';
scheme = load(fullfile(projectfolder, 'Schemes', [schemename '.mat'])).scheme;
nscheme = length(scheme);

fittingtechnique = 'LSQ';

parameter = 'fIC';
    
% Intercept in model
switch parameter
    case {'fIC', 'C'}
        intercept = false;
    case {'ADC','D', 'R', 'dEES', 'dIC', 'd', 'K', 'FA'}
        intercept = true;
end


composition = [];
paramvals = [];

for sindx = 1:length(SampleNames)

    samplename = SampleNames{sindx};
    

    % Parameter map
    paramfolder = fullfile(projectfolder, "Outputs", "Model Fitting", samplename, modeltype, schemename, fittingtechnique);
    parammap = load(fullfile(paramfolder,[parameter '.mat'])).(parameter);
    szmap = size(parammap);
    parammap_flat = reshape(parammap, [prod(szmap), 1]);
    
    % Load masks
    seriesdescription = '3DMGE_20u';
    maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, seriesdescription);
    GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
    STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
    LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;
    szmask = size(GLANDULAR);


    % == Make sample mask (COPIED FROM measure_signals.m)
    
    szbase = szmask;
    switch samplename
    
        case '20250224_UQ4'
    
            % Cylinder centred at (128, 114), radius 70 (1.4mm)
            
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-128).^2 + (Ys-114).^2 <75^2, 1, 1, szbase(3));

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;
    
        case '20250407_UQ5'
    
            % Cylinder centred at (123, 119), radius 70 (1.4mm)
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-123).^2 + (Ys-119).^2 <75^2, 1, 1, szbase(3));

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,600:end)=false;
    
    
    
        case '20250414_UQ6'
    
            % Cylinder centred at (122, 117), 
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-122).^2 + (Ys-117).^2 <75^2, 1, 1, szbase(3));

            % Remove ends
            samplemask(:,:,1:20)=false;
            samplemask(:,:,end-10:end)=false;
    
    
        case '20250522_UQ7'
    
            % Cylinder centred at (124, 117)
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-124).^2 + (Ys-117).^2 <75^2, 1, 1, szbase(3));

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;


            % EXCLUDE REGION DUE TO BUBBLES
            rows = 1:90;
            cols = 46:198;
            slices = 225:510;
            samplemask(rows, cols, slices)=false;


        case '20250523_UQ8'
    
            % Cylinder centred at (123, 119)
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-123).^2 + (Ys-119).^2 <75^2, 1, 1, szbase(3));

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;            


        case '20250524_UQ9'
    
            % Cylinder centred at (127, 118)
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-127).^2 + (Ys-118).^2 <75^2, 1, 1, szbase(3));

            % REMOVE TOP AND BOTTOM REGIONS OF MEDIUM
            samplemask(:,:,1:40)=false;
            samplemask(:,:,560:end)=false;   

    end


    % Construct composition image
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
    thiscomposition = reshape(COMPOSITION,[prod(szmap), 3]);
    bool =  sum(thiscomposition,2)>0;
    thiscomposition = thiscomposition(bool, :);
    thisparamvals = parammap_flat(bool);

    composition = [composition; thiscomposition];
    paramvals = [paramvals; thisparamvals];

end



%% Scatter plots

% x=composition(:,2);
% y=paramvals;
% 
% figure
% scatter(x, y)
% corr(x, y)
% 
% hold on;
% 
% % Linear regression line
% p = polyfit(x, y, 1);             % Fit line: y = p(1)*x + p(2)
% x_fit = linspace(min(x), max(x), 100);
% y_fit = polyval(p, x_fit);
% 
% plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Add regression line
% hold off;
% 
% xlabel('x');
% ylabel('y');
% title('Scatter Plot with Regression Line');
% grid on;

%% Linear model

% Predictors and labels
X = composition;
y = paramvals;

% Fit a multiple linear regression model
mdl = fitlm(X, y,'Intercept', intercept);

y_pred = predict(mdl, X); % Get predicted values
R2 = mdl.Rsquared.Ordinary;

Xcol = X;
Xcol(:,1)=X(:,2);
Xcol(:,2)=X(:,1);
figure;
scatter(y, y_pred,'*', MarkerEdgeAlpha=0.4, CData=Xcol); % Plot actual vs predicted
xlabel(['Actual ' parameter ' Values']);
ylabel(['Predicted ' parameter ' Values']);
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
    case {'RDI - 2 compartment - 4 param (S0)'}
        Nparam = 5;
    case {'RDI - 2 compartment - 3 param (S0)'}
        Nparam = 4;
    case {'RDI - 1 compartment - 2 param', 'DKI'}
        Nparam = 3;
    case {'ADC'}
        Nparam = 2;
end

% AIC
AIC = (nscheme)*log((resnorm)/nscheme)-2*Nparam;
AICflat = reshape(AIC, [prod(szmap), 1]);
AICvals = AICflat(bool);

figure
histogram(AICvals)


 %% Non-linear models

% 
% % X: input matrix (nFeatures × nSamples)
% % y: output vector (1 × nSamples)
% hiddenLayerSize = 10;
% net = fitnet(hiddenLayerSize);  % e.g., 
% net = train(net, transpose(X), transpose(y));
% 
% 
% y_pred = transpose(net(transpose(X))); 
% 
% figure
% scatter(y, y_pred)