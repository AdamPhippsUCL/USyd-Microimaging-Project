% Script to measure S, G, and L signals

clear;
projectfolder = pwd;


%% Sample and image details

% Sample
samplename = '20250224_UQ4';

% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end


% == Base image (3DMGE)

seriesdescription = '3DMGE_20u';
baseImg = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', samplename, seriesdescription, 'avgImageArray.mat')).avgImageArray;
szbase = size(baseImg);


% == Low resolution DW-MRI images

schemename = '20250224_UQ4 AllDELTA';
schemesfolder = fullfile(projectfolder, 'Schemes');
load(fullfile(schemesfolder, schemename));
nscheme = length(scheme);

% Correct for non-zero b0 value
effb0 = 10.3;
for indx = 1:nscheme
   if scheme(indx).bval==0
       continue
   end
   scheme(indx).bval = scheme(indx).bval-effb0;
end

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

Nimg = nscheme;


%% Load masks

maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, seriesdescription);
GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;


% displaymasks = zeros([size(GLANDULAR), 3]);
% displaymasks(:,:,:,1) = logical(GLANDULAR);
% displaymasks(:,:,:,2) = logical(STROMA);
% displaymasks(:,:,:,3) = logical(LUMEN);
% 
% sl = 128;
% figure
% imshow(squeeze(baseImg(sl,:,:)),[]);
% hold on
% mask = imshow(squeeze(displaymasks(sl,:,:,:)));
% set(mask, 'AlphaData', 0.2)


%% Diffuions imaging data preprocessing

ImageArrays = struct();
DINFOS = struct();

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    thisfolder = fullfile(ImagingDataFolder, samplename, SeriesDescription);

    ImageArray = load(fullfile(thisfolder, 'axialImageArray.mat')).ImageArray;
    dinfo = load(fullfile(thisfolder, 'axialdinfo.mat')).dinfo;

    % Append to structures
    DINFOS(seriesindx).dinfo = dinfo;
    ImageArrays(seriesindx).ImageArray = ImageArray;

end

szimg = size(ImageArray, 1:3);
IMGS = ones([Nimg, szimg]);

% seriesindx = 1 used for data normalisation
img = ImageArrays(1).ImageArray;
dinfo = DINFOS(1).dinfo;
bvals = [dinfo(:).DiffusionBValue];
b0bools = (bvals==0);
b0imgs = img(:,:,:,b0bools);
b0img = mean(b0imgs,4);
meanb0 = mean(b0img(:));


% seriesindx > 1 used for diffusion data
for seriesindx = 2:length(SeriesDescriptions)

    % Load image
    img = ImageArrays(seriesindx).ImageArray;

    % Load dinfo
    dinfo = DINFOS(seriesindx).dinfo;
    bvals = [dinfo(:).DiffusionBValue];

    % b0 imgs
    b0bools = (bvals==0);
    b0imgs = img(:,:,:,b0bools);
    thisb0img = mean(b0imgs,4);
    thismeanb0 = mean(thisb0img(:));

    % b imgs
    bbools = (bvals > 0);
    bimgs = img(:,:,:,bbools);
    bimg = mean(bimgs,4);
    bval = bvals(find(bbools,1));

    % Normalize and append to Y array
    IMGS(seriesindx,:,:,:) = (meanb0/thismeanb0)*(bimg./b0img); 

end

clear DINFOS
clear ImageArrays

%% Make sample mask

switch samplename

    case '20250224_UQ4'

        % Cylinder centred at (128, 114)
        
        samplemask = zeros(size(baseImg));
        
        [Xs, Ys] = meshgrid( ...
           1:size(baseImg,2), ...
           1:size(baseImg,1) ...
            );

        samplemask(:,:,:) = repmat((Xs-128).^2 + (Ys-114).^2 <75^2, 1, 1, size(baseImg,3));

end




%% Construct composition map

% Map size
szmap = szimg;
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


% Select voxels with non-zero composition

composition = reshape(COMPOSITION,[prod(szmap), 3]);
bool =  sum(composition,2)>0;
composition = composition(bool, :);

imgs = reshape(IMGS, [Nimg, prod(szimg)]);
imgs = imgs(:,bool);



%% SGL signal measurement (Linear regression)

% Predictors [fs, fg, fl]
X = composition;

% Linear regression function
func = @(b, X) signal_func(b, X);

% Initial guess and bounds [Ss, Sg, Sl]
beta0 = [0.5, 0.5, 0.5];
lb = [0,0,0];
ub=[1,1,1];

% Lumen diffusivity
Dl = 0.002;
err = 0.0001;
% Initialise array for signal measurements
signals = zeros(3, Nimg, 2);

for imgindx = 1:Nimg

    if imgindx == 1
        signals(:,imgindx,1)=1;
        signals(:,imgindx,2)=0;
        continue
    end

    y = transpose(imgs(imgindx, :));

    % Set lumen signal
    Sl = exp(-scheme(imgindx).bval*Dl);
    lb(3)=Sl-err;
    ub(3)=Sl+err;

    % Apply fitting
    [beta_fit, resnorm, residual, exitflag, output, lambda, jacobian]=lsqcurvefit(func, beta0, X, y, lb, ub);

    % Standard error calculation
    mse = resnorm / (length(y) - length(beta0));
    cov_matrix = mse * inv(jacobian' * jacobian);
    stderr = sqrt(diag(cov_matrix));  


    signals(:,imgindx,1) = beta_fit;
    signals(:,imgindx,2) = stderr;


end



%% Display measured signals

figure
bshift = 10;

% Short Delta

indices = 2:6;
s = signals(:,indices,:);
bvals = [scheme(indices).bval];

errorbar(bvals-1*bshift, s(1,:,1), s(1,:,2), '--*', color=[0.4660 0.6740 0.1880], DisplayName='Stroma (Short Delta)');
hold on
errorbar(bvals-1*bshift, s(2,:,1), s(2,:,2), '--*', color=[0.8500 0.3250 0.0980], DisplayName='Glandular (Short Delta)');
errorbar(bvals-1*bshift, s(3,:,1), s(3,:,2), '--*', color=[0 0.4470 0.7410], DisplayName = 'Lumen (Short Delta)');

% Long Delta
indices = 7:11;
s = signals(:,indices,:);

errorbar(bvals+1*bshift, s(1,:,1), s(1,:,2), '--.', color=[0.4660 0.6740 0.1880], DisplayName='Stroma (Long Delta)');
errorbar(bvals+1*bshift, s(2,:,1), s(2,:,2), '--.', color=[0.8500 0.3250 0.0980], DisplayName='Glandular (Long Delta)');
errorbar(bvals+1*bshift, s(3,:,1), s(3,:,2), '--.', color=[0 0.4470 0.7410], DisplayName = 'Lumen (Long Delta)');

xticks(bvals); 
xticklabels(bvals)
ylim([-0.1,0.8])
xlim([800,2200])
xlabel('b-value')
xticks(bvals)
xticklabels(["1000", "1250", "1500", "1750", "2000"])
grid on
legend;
ax=gca();
ax.FontSize=14;


%% Save measured signals

savefolder = fullfile(projectfolder, "Outputs", "Signal Measurement", samplename, schemename);
mkdir(savefolder);
save(fullfile(savefolder, 'signals.mat'), 'signals')
save(fullfile(savefolder, 'scheme.mat'), 'scheme')

%% Function

function signals = signal_func(b,X)

% b = [Ss, Sg, Sl]
% X = [fs1, fg1, fl1; fs2, fg2, fl2; ... ]

Ss = b(1);
Sg = b(2);
Sl = b(3);

N = size(X,1);
signals = zeros(N,1);

for indx = 1:N

    fs = X(indx,1);
    fg = X(indx,2);
    fl = X(indx,3);

    signals(indx) = fs*Ss + fg*Sg + fl*Sl;
    
end


end