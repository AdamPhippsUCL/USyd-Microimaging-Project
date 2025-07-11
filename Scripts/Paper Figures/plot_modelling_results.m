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


ShowErrorMap = false;

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
            fit_AIC = load(fullfile(outputfolder, 'AIC.mat')).AIC;

            % Predicted model parameters
            pred_fIC = S_params(1).*sample_composition(:,:,:,1) + E_params(1).*sample_composition(:,:,:,2);
            pred_R = S_params(2).*sample_composition(:,:,:,1) + E_params(2).*sample_composition(:,:,:,2) + 6.5.*sample_composition(:,:,:,3);
            pred_dIC = S_params(3).*sample_composition(:,:,:,1) + E_params(3).*sample_composition(:,:,:,2) + L_params(2).*sample_composition(:,:,:,3);
            pred_dEES = S_params(4).*sample_composition(:,:,:,1) + E_params(4).*sample_composition(:,:,:,2) + L_params(2).*sample_composition(:,:,:,3);

    end

    szmap = size(pred_fIC);

    this_composition = reshape(sample_composition,[prod(size(sample_composition, 1:3)), 3]);
    bool =  sum(this_composition,2)>0;
    this_composition = this_composition(bool, :);
    bool = reshape(bool, size(sample_composition, 1:3));

    % Inspect differences in pred and fit maps here...
    if ShowErrorMap
        % Plot Estimated - predicted error over MGE as in other script
        MGE_seriesdescription = '3DMGE_20u';
        MGE = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, MGE_seriesdescription, 'avgImageArray.mat')).avgImageArray;
    
        fIC_error = (fit_fIC - pred_fIC).*bool;
        fIC_error_large = zeros(size(MGE));
    
        ResFactor = size(MGE)./size(fIC_error);
    
        for rindx = 1:szmap(1)
            for cindx = 1:szmap(2)
                for slindx = 1:szmap(3)
                    rows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
                    cols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
                    slices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));
                    fIC_error_large(rows,cols,slices)=fIC_error(rindx, cindx, slindx);  
                end
            end
        end
    
        [Xs, Ys] = meshgrid( ...
            linspace(1,szmap(3),size(MGE,3)), ...
            linspace(1,szmap(2),size(MGE,2)) ...
            );    
    
    
        xs = 1:240;
        ys = 1:640;
        sl = 120;
        f=figure;
    
        ax1 = axes;
        pcolor(ax1, Xs(xs, ys), Ys(xs, ys), rescale(squeeze(MGE(sl,xs,ys)), 0,1))
        shading flat; % Remove grid-like shading
        grid off;
        colormap(ax1,gray);
        xticks([])
        yticks([])
        current_pos = ax1.Position;
        ax1.Position = current_pos - [0.1, 0, 0, 0];
        hold on
    
        ax2 = axes;
        m=pcolor(ax2, Xs(xs, ys), Ys(xs, ys), squeeze(fIC_error_large(sl,xs,ys)) );
        caxis([-0.16, 0.16])
        % caxis([-max(abs(error(:))) max(abs(error(:)))]);
        shading flat; % Remove grid-like shading
        grid off;
        set(m, 'FaceAlpha', 0.3);
        linkaxes([ax1, ax2]);
        ax2.Visible = 'off'; % Hide second axes to avoid overlapping ticks
        colormap(ax2,redblue);
        c=colorbar;
        c.Label.String = 'Sphere fraction error';
        ax2.Position = ax1.Position;
    
        set(ax1, 'YDir', 'reverse');
        set(ax2, 'YDir', 'reverse'); 
        % f.Position = [0.0010    0.0490    1.5360    0.7408]*1e3;
        ax1.FontSize = 12;
        ax2.FontSize = 14;
    end
    

    % AIC for non-zero composition voxels
    this_fit_AIC = fit_AIC(bool);

    switch ModelName

        case 'RDI - 2 compartment - 4 param (S0)'

        % Extract modelling and prediction results for non-zero composition
        this_fit_fIC = fit_fIC(bool);
        this_fit_R = fit_R(bool);
        this_fit_dIC = fit_dIC(bool);
        this_fit_dEES = fit_dEES(bool);

        % disp(max(this_fit_fIC))
    
        this_pred_fIC = pred_fIC(bool);
        this_pred_R = pred_R(bool);
        this_pred_dIC = pred_dIC(bool);
        this_pred_dEES = pred_dEES(bool);
    
        % Append to parameter array
        if sampleindx == 1       
            composition = this_composition;
            pred_params = [this_pred_fIC'; this_pred_R'; this_pred_dIC'; this_pred_dEES'];
            fit_params = [this_fit_fIC'; this_fit_R'; this_fit_dIC'; this_fit_dEES'; this_fit_AIC'];
        else
            composition = cat(1, composition, this_composition);
            pred_params = cat(2, pred_params, [this_pred_fIC'; this_pred_R'; this_pred_dIC'; this_pred_dEES']);
            fit_params = cat(2, fit_params, [this_fit_fIC'; this_fit_R'; this_fit_dIC'; this_fit_dEES'; this_fit_AIC']);
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
scatter(pred_params(1,:), fit_params(1,:),  '*', MarkerEdgeAlpha=0.4, CData=new_composition)
hold on
plot([0, 0.3], [0, 0.3], 'k')
grid on
xlim([-0.025,0.325])
ylim([-0.05,0.6])
xlabel('Predicted sphere fraction')
ylabel('Estimated sphere fraction')

% R2 value
SSres = sum( (pred_params(1,:) - fit_params(1,:)).^2 );
SStot = length(fit_params(1,:))*var(fit_params(1,:));
R2 = 1-SSres/SStot;

text(0.03, 0.945, ['R^2 = ' sprintf( '%0.3f', R2) ], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border


% figure
% scatter(pred_params(2,:), fit_params(2,:),  '*', MarkerEdgeAlpha=0.4, CData=new_composition)
% grid on
% xlim([0,10])
% ylim([0,10])
% 
% figure
% scatter(pred_params(3,:), fit_params(3,:),  '*', MarkerEdgeAlpha=0.4, CData=new_composition)
% grid on
% xlim([0,2])
% ylim([0,2])


figure
scatter(pred_params(4,:), fit_params(4,:),  '*', MarkerEdgeAlpha=0.4, CData=new_composition)
hold on
plot([0.5, 2], [0.5, 2], 'k')
grid on
xlim([0.35,2.05])
ylim([0.35,2.05])
xlabel('Predicted d_{out} (s/mm^2)')
ylabel('Estimated d_{out} (s/mm^2)')

% R2 value
SSres = sum( (pred_params(4,:) - fit_params(4,:)).^2 );
SStot = length(fit_params(4,:))*var(fit_params(4,:));
R2 = 1-SSres/SStot;

text(0.03, 0.945, ['R^2 = ' sprintf( '%0.3f', R2) ], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border




% AIC
figure
scatter(fit_params(1,:), fit_params(5,:),  '*', MarkerEdgeAlpha=0.4, CData=new_composition)
grid on
xlim([-0.02, 0.57])
ylim([-175, -95])
xlabel('Estimated sphere fraction')
ylabel('AIC')

%%


function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b]; 


end