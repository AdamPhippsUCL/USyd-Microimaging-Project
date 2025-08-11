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

ModelName = 'Ball+Sphere';
schemename = '20250224_UQ4 AllDELTA';
fittingtechnique = 'LSQ';

% Load ESL modelling estimates
ESL_Model_RESULTS = load(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'Modelling', 'RESULTS.mat')).RESULTS;
S_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'S'), strcmp({ESL_Model_RESULTS(:).ModelName}, ModelName))).ModelParams;
E_params = ESL_Model_RESULTS(and(strcmp({ESL_Model_RESULTS(:).Component}, 'E'), strcmp({ESL_Model_RESULTS(:).ModelName}, ModelName))).ModelParams;
L_params = [1, 2]; % ADC parameters...

ShowErrorMap = false;

% Load model parameters from direct fitting, and composition array
for sampleindx = 1:length(SampleNames)

    SampleName = SampleNames{sampleindx};

    % Load composition array
    sample_composition = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;
  
    % Cancer samples = UQ4B, UQ4M, UQ6N
    switch SampleName(end-2:end)
        case 'UQ4'
            % Only include sample N
            NMASK = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'NMASK.mat')).NMASK;
            sample_composition = sample_composition.*double(NMASK);
            clear NMASK
        case 'UQ6'
            % Only include samples B and M
            BMASK = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'BMASK.mat')).BMASK;
            MMASK = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'MMASK.mat')).MMASK;
            sample_composition = sample_composition.*(double(BMASK)+double(MMASK));
    end



    switch ModelName
        
        case 'Ball+Sphere'

            outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, ModelName, schemename, fittingtechnique);
            fit_fs = load(fullfile(outputfolder, 'fs.mat')).fs;
            fit_Ds = load(fullfile(outputfolder, 'Ds.mat')).Ds;
            fit_R = load(fullfile(outputfolder, 'R.mat')).R;
            fit_Db = load(fullfile(outputfolder, 'Db.mat')).Db;
            fit_AIC = load(fullfile(outputfolder, 'AIC.mat')).AIC;

            % Predicted model parameters
            pred_fs = S_params(1).*sample_composition(:,:,:,2) + E_params(1).*sample_composition(:,:,:,1);
            pred_R = S_params(2).*sample_composition(:,:,:,2) + E_params(2).*sample_composition(:,:,:,1) + 6.5.*sample_composition(:,:,:,3);
            pred_Ds = S_params(3).*sample_composition(:,:,:,2) + E_params(3).*sample_composition(:,:,:,1) + L_params(2).*sample_composition(:,:,:,3);
            pred_Db = S_params(4).*sample_composition(:,:,:,2) + E_params(4).*sample_composition(:,:,:,1) + L_params(2).*sample_composition(:,:,:,3);

    end

    szmap = size(pred_fs);

    this_composition = reshape(sample_composition,[prod(size(sample_composition, 1:3)), 3]);
    bool =  sum(this_composition,2)>0;
    this_composition = this_composition(bool, :);
    bool = reshape(bool, size(sample_composition, 1:3));

    % Inspect differences in pred and fit maps here...
    if ShowErrorMap
        % Plot Estimated - predicted error over MGE as in other script
        MGE_seriesdescription = '3DMGE_20u';
        MGE = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, MGE_seriesdescription, 'avgImageArray.mat')).avgImageArray;
    
        fs_error = (fit_fs - pred_fs).*bool;
        fs_error_large = zeros(size(MGE));
    
        ResFactor = size(MGE)./size(fs_error);
    
        for rindx = 1:szmap(1)
            for cindx = 1:szmap(2)
                for slindx = 1:szmap(3)
                    rows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
                    cols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
                    slices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));
                    fs_error_large(rows,cols,slices)=fs_error(rindx, cindx, slindx);  
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
        m=pcolor(ax2, Xs(xs, ys), Ys(xs, ys), squeeze(fs_error_large(sl,xs,ys)) );
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

        case 'Ball+Sphere'

        % Extract modelling and prediction results for non-zero composition
        this_fit_fs = fit_fs(bool);
        this_fit_R = fit_R(bool);
        this_fit_Ds = fit_Ds(bool);
        this_fit_Db = fit_Db(bool);
    
        this_pred_fs = pred_fs(bool);
        this_pred_R = pred_R(bool);
        this_pred_Ds = pred_Ds(bool);
        this_pred_Db = pred_Db(bool);
    
        % Append to parameter array
        if sampleindx == 1       
            composition = this_composition;
            pred_params = [this_pred_fs'; this_pred_R'; this_pred_Ds'; this_pred_Db'];
            fit_params = [this_fit_fs'; this_fit_R'; this_fit_Ds'; this_fit_Db'; this_fit_AIC'];
        else
            composition = cat(1, composition, this_composition);
            pred_params = cat(2, pred_params, [this_pred_fs'; this_pred_R'; this_pred_Ds'; this_pred_Db']);
            fit_params = cat(2, fit_params, [this_fit_fs'; this_fit_R'; this_fit_Ds'; this_fit_Db'; this_fit_AIC']);
        end
        
    end


end

%% Display results

f1=figure;
scatter(pred_params(1,:), fit_params(1,:),  6, 'filled', 'MarkerFaceAlpha', 0.7, CData= composition);
hold on
plot([0, 0.3], [0, 0.3], 'k')
grid on
xlim([-0.025,0.325])
ylim([-0.05,0.48])
xlabel('Predicted Sphere Fraction')
ylabel('Estimated Sphere Fraction')

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

ax = gca();
ax.FontSize = 12;

saveas(f1, fullfile(projectfolder, 'Figures', 'Predicted vs Estimated Sphere Fraction.png'))




f2=figure;
scatter(pred_params(4,:), fit_params(4,:), 6, 'filled', 'MarkerFaceAlpha', 0.7, CData= composition);
hold on
plot([0.5, 2], [0.5, 2], 'k')
grid on
xlim([0.35,2.2])
ylim([0.35,2.2])
xlabel('Predicted D_{b} (x10^{-3} mm^2/s)')
ylabel('Estimated D_{b} (x10^{-3} mm^2/s)')

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

ax = gca();
ax.FontSize = 12;

saveas(f2, fullfile(projectfolder, 'Figures', 'Predicted vs Estimated Db.png'))

% AIC
figure
scatter(fit_params(1,:), fit_params(5,:), 6, 'filled', 'MarkerFaceAlpha', 0.7, CData= composition);
grid on
xlim([-0.02, 0.57])
ylim([-175, -95])
xlabel('Estimated sphere fraction')
ylabel('AIC')


% BLAND ALTMAN

fs_avg = (pred_params(1,:)+fit_params(1,:))/2;
fs_diff = fit_params(1,:)-pred_params(1,:);
fs_LOA = [mean(fs_diff), mean(fs_diff)-1.96*std(fs_diff), mean(fs_diff)+1.96*std(fs_diff)];
save(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'Modelling', 'fs_LOA.mat'), 'fs_LOA')

f3 = figure;
scatter(fs_avg, fs_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=  composition, HandleVisibility='off');
yline(mean(fs_diff), '-', DisplayName='Bias', LineWidth=1)
hold on
yline(mean(fs_diff)+1.96*std(fs_diff), '-.', DisplayName='95% LOA', LineWidth=1)
yline(mean(fs_diff)-1.96*std(fs_diff), '-.', HandleVisibility="off", LineWidth=1)
xlim([-0.02 0.37])
ylim([-0.35, 0.37])
xlabel('Mean of Predicted and Estimated Sphere Fraction')
ylabel('Estimated - Predicted Sphere Fraction')
legend
grid on
ax = gca();
ax.FontSize = 12;
saveas(f3, fullfile(projectfolder, 'Figures', ['Bland Altman Benign Sphere Fraction.png']))


Db_avg = (pred_params(4,:)+fit_params(4,:))/2;
Db_diff = fit_params(4,:)-pred_params(4,:);
Db_LOA = [mean(Db_diff), mean(Db_diff)-1.96*std(Db_diff), mean(Db_diff)+1.96*std(Db_diff)];
save(fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation', 'Multi-sample', 'Modelling', 'Db_LOA.mat'), 'Db_LOA')

f4 = figure;
scatter(Db_avg, Db_diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=  composition, HandleVisibility='off');
yline(mean(Db_diff), '-', DisplayName='Bias', LineWidth=1)
hold on
yline(mean(Db_diff)+1.96*std(Db_diff), '-.', DisplayName='95% LOA', LineWidth=1)
yline(mean(Db_diff)-1.96*std(Db_diff), '-.', HandleVisibility="off", LineWidth=1)
xlim([0.4 2.12])
ylim([-1.1, 1.1])
xlabel('Mean of Predicted and Estimated D_{b} (x10^{-3} mm^2/s)')
ylabel('Estimated - Predicted D_{b} (x10^{-3} mm^2/s)')
legend
grid on
ax = gca();
ax.FontSize = 12;
saveas(f4, fullfile(projectfolder, 'Figures', ['Bland Altman Benign Dout.png']))


% %%
% 
% function c = redblue(m)
% %REDBLUE    Shades of red and blue color map
% %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
% %   The colors begin with bright blue, range through shades of
% %   blue to white, and then through shades of red to bright red.
% %   REDBLUE, by itself, is the same length as the current figure's
% %   colormap. If no figure exists, MATLAB creates one.
% %
% %   For example, to reset the colormap of the current figure:
% %
% %             colormap(redblue)
% %
% %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
% %   COLORMAP, RGBPLOT.
% %   Adam Auton, 9th October 2009
% if nargin < 1, m = size(get(gcf,'colormap'),1); end
% if (mod(m,2) == 0)
%     % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
%     m1 = m*0.5;
%     r = (0:m1-1)'/max(m1-1,1);
%     g = r;
%     r = [r; ones(m1,1)];
%     g = [g; flipud(g)];
%     b = flipud(r);
% else
%     % From [0 0 1] to [1 1 1] to [1 0 0];
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
% end
% c = [r g b]; 


end