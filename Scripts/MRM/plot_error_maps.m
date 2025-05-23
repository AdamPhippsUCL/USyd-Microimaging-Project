% Script to visualise error of linear model across sample

clear;
projectfolder = pwd;

%% Sample and image details

% Sample
SampleName = '20250224_UQ4'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6'


% Image
seriesindx = 5;
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

SeriesDescription = SeriesDescriptions{seriesindx};

ImageFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, SeriesDescription);
ImageArray = load(fullfile(ImageFolder,'normalisedImageArray.mat')).ImageArray;
% dinfo = load(fullfile(ImageFolder,'axialdinfo.mat')).dinfo;




%% Load signal measurements and composition fractions

multisample = true;

switch multisample

    case true
        signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'signals.mat')).signals;
        signals = squeeze(signals(:,seriesindx,1));

    case false
        signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', SampleName, 'signals.mat')).signals;
        signals = squeeze(signals(:,seriesindx,1));
end

COMPOSITION = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;
        


%% Calculate predicted normalised image and error

signals = reshape(signals, [1,1,1,3]);
pred = sum(COMPOSITION.*repmat(signals, [size(COMPOSITION, 1:3)]), 4);
error = (pred - ImageArray).*double(pred>0);
szmap = size(pred);


%% Load high-res MGE image

MGE_seriesdescription = '3DMGE_20u';
MGE = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, MGE_seriesdescription, 'avgImageArray.mat')).avgImageArray;

%% Display error on top of MGE

ResFactor = size(MGE)./size(error);

error_large = zeros(size(MGE));

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
        for slindx = 1:szmap(3)

            rows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
            cols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
            slices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));

            error_large(rows,cols,slices)=error(rindx, cindx, slindx);

        end
    end
end


[Xs, Ys] = meshgrid( ...
    linspace(1,szmap(3),size(MGE,3)), ...
    linspace(1,szmap(2),size(MGE,2)) ...
    );

sl = 120;

figure
ax1 = axes;
pcolor(ax1, Xs, Ys, rescale(squeeze(MGE(sl,:,:)), 0,1))
shading flat; % Remove grid-like shading
grid off;
colormap(ax1,gray);


hold on
ax2 = axes;
m=pcolor(ax2, Xs, Ys, squeeze(error_large(sl,:,:)) );
caxis([-0.15, 0.15])
% caxis([-max(abs(error(:))) max(abs(error(:)))]);
shading flat; % Remove grid-like shading
grid off;
set(m, 'FaceAlpha', 0.3);
linkaxes([ax1, ax2]);
ax2.Visible = 'off'; % Hide second axes to avoid overlapping ticks
colormap(ax2,redblue);
colorbar;
ax2.Position = ax1.Position;
% 

set(ax1, 'YDir', 'reverse');
set(ax2, 'YDir', 'reverse');











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