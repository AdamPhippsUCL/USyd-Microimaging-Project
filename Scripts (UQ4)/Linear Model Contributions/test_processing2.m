
clear;

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

samplename = '20250224_UQ4';

%% Base image (3DMGE)

seriesdescription = '3DMGE_20u';

% Load image
baseImg = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', samplename, seriesdescription, 'avgImageArray.mat')).avgImageArray;
szbase = size(baseImg);


%% MASKS

maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, seriesdescription);
GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;


%% Display masks

displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = logical(GLANDULAR);
displaymasks(:,:,:,2) = logical(STROMA);
displaymasks(:,:,:,3) = logical(LUMEN);

sl = 128;
figure
imshow(squeeze(baseImg(sl,:,:)),[]);
hold on
mask = imshow(squeeze(displaymasks(sl,:,:,:)));
set(mask, 'AlphaData', 0.2)


%% Load parameter maps

% Model type
% modeltype = 'ADC';
modeltype = 'RDI - 2 compartment - 4 param';

% Parameter
parameter = 'C';

% Scheme name
schemename = '20250224_UQ4 AllDELTA';
% schemename = 'STEAM_ShortDELTA_50 (640 micron)';

% Fitting technique
fittingtechnique = 'LSQ';


% Parameter folder
paramfolder = fullfile(projectfolder, "Outputs", "Model Fitting", samplename, modeltype, schemename, fittingtechnique);

% Load parameter map
switch parameter
    case 'ADC'
        parammap = load(fullfile(paramfolder, 'ADC.mat')).ADC;
        intercept = true;
    case 'fIC'
        parammap = load(fullfile(paramfolder, 'fIC.mat')).fIC;
        intercept = false;
    case 'D'
        parammap = load(fullfile(paramfolder, 'D.mat')).D;
        intercept = true;
    case 'R'
        parammap = load(fullfile(paramfolder, 'R.mat')).R;
        intercept = true;
    case 'C'
        parammap = load(fullfile(paramfolder, 'C.mat')).C;
        intercept = false;
    case 'dEES'
        parammap = load(fullfile(paramfolder, 'dEES.mat')).dEES;
        intercept = true;
    case 'dIC'
        parammap = load(fullfile(paramfolder, 'dIC.mat')).dIC;
        intercept = true;
    case 'd'
        parammap = load(fullfile(paramfolder, 'd.mat')).d;
        intercept = true;
    case 'K'
        parammap = load(fullfile(paramfolder, 'K.mat')).K;
        intercept = true;        
    case 'FA'
        parammap = load(fullfile(paramfolder, 'FA.mat')).FA;
    case 'fs'
        parammap = load(fullfile(paramfolder, 'fs.mat')).fs;
    case 'fg'
        parammap = load(fullfile(paramfolder, 'fg.mat')).fg;
    case 'fl'
        parammap = load(fullfile(paramfolder, 'fl.mat')).fl;
end



szmap = size(parammap);



%% Make sample mask


switch samplename

    case '20250224_UQ4'

        % Cylinder centred at (128, 114)
        
        samplemask = zeros(size(baseImg));
        
        [Xs, Ys] = meshgrid( ...
           1:size(baseImg,2), ...
           1:size(baseImg,1) ...
            );
        
        samplemask(:,:,:) = repmat((Xs-128).^2 + (Ys-114).^2 < 75^2, 1, 1, size(baseImg,3));

end
% for indx = 1:size(baseImg, 1)
%     vals = squeeze(baseImg(indx,:,300));
%     smoothvals = sgolayfilt(vals, 1, 11);
% 
%     diffs = diff(smoothvals);
%     [~,xmin] = max(diffs);
%     [~,xmax] = min(diffs);
%     xmin = xmin; %+ 10;
%     xmax = xmax; %- 10;
% 
%     samplemask(indx, xmin:xmax,:) = 1;
% 
% 
% end

% % Test slice
% sl = 140;
% figure
% imshow(squeeze(baseImg(sl,:,:)),[]);
% hold on
% m=imshow(squeeze(max(baseImg(:))*samplemask(sl,:,:)),[]);
% set(m, 'AlphaData', 0.5)


%% Construct composition image

COMPOSITION = zeros([szmap, 3]);
ResFactor = szbase./szmap;
Nvoxel = prod(ResFactor);

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
        for slindx = 1:szmap(3)

            baserows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
            basecols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
            baseslices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));

            % Test in sample
            if ~all(logical(samplemask(baserows, basecols, baseslices)), "all")
                continue
            end

            % COMPOSITION         
            COMPOSITION(rindx, cindx, slindx, 1) = sum(double(STROMA(baserows, basecols, baseslices)), "all")/Nvoxel;
            COMPOSITION(rindx, cindx, slindx, 2) = sum(double(GLANDULAR(baserows, basecols, baseslices)), "all")/Nvoxel;
            COMPOSITION(rindx, cindx, slindx, 3) = sum(double(LUMEN(baserows, basecols, baseslices)), "all")/Nvoxel;        

        end
    end
end

STROMAvals = squeeze(COMPOSITION(:,:,:,1));
GLANDULARvals = squeeze(COMPOSITION(:,:,:,2));
LUMENvals = squeeze(COMPOSITION(:,:,:,3));

%% TEST SCATTER


% Param vals
paramvals = parammap(:);

stromavals = STROMAvals(:);
glandularvals = GLANDULARvals(:);
lumenvals = LUMENvals(:);

mask = (stromavals + glandularvals + lumenvals > 0);

figure
scatter(glandularvals(mask),parammap(mask) )
lsline
title('G')
% 
% figure
% scatter(stromavals(mask),parammap(mask) )
% lsline
% title('S')
% 
% figure
% scatter(lumenvals(mask),parammap(mask) )
% lsline
% title('L')


%% Multiple linear model

X = zeros([length(parammap(mask)), 3]);
X(:,1) = stromavals(mask);
X(:,2) = glandularvals(mask);
X(:,3) = lumenvals(mask);

Y = paramvals(mask);



% Fit a multiple linear regression model
mdl = fitlm(X, Y,'Intercept', intercept);

% Display regression summary
disp(mdl);


Y_pred = predict(mdl, X); % Get predicted values
r = corr(Y, Y_pred);

figure;
scatter(Y, Y_pred, 'filled'); % Plot actual vs predicted
xlabel('Actual Values');
ylabel('Predicted Values');
title('Actual vs. Predicted Values (LM)');
hold on;
plot(Y, Y, 'r--'); % 45-degree reference line

switch intercept
    case true
        b = mdl.Coefficients.Estimate(2:end); % Exclude intercept

        sx = std(X);  % Standard deviation of predictors
        sy = std(Y);  % Standard deviation of response
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
        sy = std(Y);  % Standard deviation of response
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

% %% Non-linear model
% 
% nlmodel = @(b,X) b(1)*X(:,1) + b(2)*(X(:,1).^2) + b(3)*X(:,2) + b(4)*(X(:,2).^2) + b(5)*X(:,3) + b(6)*(X(:,3).^2);
% 
% beta0 = [0,0,0,0,0,0];
% 
% mdl = fitnlm(X, Y, nlmodel, beta0);
% disp(mdl);
% 
% Y_pred = predict(mdl, X); % Get predicted values
% r = corr(Y, Y_pred);
% 
% figure;
% scatter(Y, Y_pred, 'filled'); % Plot actual vs predicted
% xlabel('Actual Values');
% ylabel('Predicted Values');
% title('Actual vs. Predicted Values (NLM)');
% hold on;
% plot(Y, Y, 'r--'); % 45-degree reference line


%% Display parameter map as overlay

param_large = zeros(size(baseImg));

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
        for slindx = 1:szmap(3)

            baserows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
            basecols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
            baseslices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));

            % Test in sample
            if ~all(logical(samplemask(baserows, basecols, baseslices)), "all")
                continue
            end

            param_large(baserows,basecols,baseslices)=parammap(rindx, cindx, slindx);

        end
    end
end


[Xsmall, Ysmall] = meshgrid( ...
    linspace(1,size(parammap,3),size(baseImg,3)), ...
    linspace(1,size(parammap,2),size(baseImg,2)) ...
    );

figure
ax1 = axes;
pcolor(ax1, Xsmall, Ysmall, rescale(squeeze(baseImg(120,:,:)), 0,1))
shading flat; % Remove grid-like shading
grid off;
colormap(ax1,gray);


hold on
ax2 = axes;
m=pcolor(ax2, Xsmall, Ysmall, rescale(squeeze(param_large(120,:,:)), 0,1));
shading flat; % Remove grid-like shading
grid off;
set(m, 'FaceAlpha', 0.5);
linkaxes([ax1, ax2]);
ax2.Visible = 'off'; % Hide second axes to avoid overlapping ticks
colormap(ax2,hot);
% 

set(ax1, 'YDir', 'reverse');
set(ax2, 'YDir', 'reverse');

% % hold on
% % ax2 = axes;
% % m=pcolor(ax2,Xlarge, Ylarge, rescale(squeeze(parammap(6,:,:)), 0,1));
% 
% 
% shading flat; % Remove grid-like shading
% % colorbar;
% grid off;
% set(m, 'FaceAlpha', 0.5);
% colormap(ax2,jet);
% daspect([2 5 1])
% linkaxes([ax1, ax2]);
% ax2.Visible = 'off'; % Hide second axes to avoid overlapping ticks
% 
% 
% 
