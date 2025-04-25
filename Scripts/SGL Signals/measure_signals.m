% Script to measure S, G, and L signals

clear;
projectfolder = pwd;


%% Sample and image details

% Sample
SampleNames = {'20250224_UQ4', '20250407_UQ5'};%, '20250414_UQ6'};
% SampleNames = {'20250414_UQ6'};

schemename = '20250224_UQ4 AllDELTA';
schemesfolder = fullfile(projectfolder, 'Schemes');
load(fullfile(schemesfolder, schemename));
nscheme = length(scheme);
Nimg = nscheme;
effb0vals = [scheme(:).effb0vals];

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


% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end
    

composition = [];
imgs = [];

for sindx = 1:length(SampleNames)

    samplename = SampleNames{sindx};

    % == Load masks

    baseseriesdescription = '3DMGE_20u';
    maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, baseseriesdescription);
    GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
    STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
    LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;
    szbase = size(GLANDULAR);


    % == Diffuion imaging data preprocessing
    
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
    
    b0prctiles = [
        prctile(b0img(:), 10),...
        prctile(b0img(:), 20),...
        prctile(b0img(:), 30),...
        prctile(b0img(:), 40),...
        prctile(b0img(:), 50),...
        prctile(b0img(:), 60),...
        prctile(b0img(:), 70),...
        prctile(b0img(:), 80),...
        prctile(b0img(:), 90),...
    ];
   
    
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
    
        thisb0prctiles = [
            prctile(thisb0img(:), 10),...
            prctile(thisb0img(:), 20),...
            prctile(thisb0img(:), 30),...
            prctile(thisb0img(:), 40),...
            prctile(thisb0img(:), 50),...
            prctile(thisb0img(:), 60),...
            prctile(thisb0img(:), 70),...
            prctile(thisb0img(:), 80),...
            prctile(thisb0img(:), 90),...
        ];

        p = polyfit(b0prctiles, thisb0prctiles, 1);  


        % b imgs
        bbools = (bvals > 0);
        bimgs = img(:,:,:,bbools);
        bimg = mean(bimgs,4);
        bval = bvals(find(bbools,1));

        normimg = (1/p(1))*(bimg./b0img);
        IMGS(seriesindx,:,:,:) =  normimg;


        % % Correct for eff b0 difference using high res D...
        % D=load(fullfile(projectfolder, 'Scripts', 'SGL Signals', 'D', 'SampleName', 'D_DS')).D_DS;
        % effb0val = effb0vals(seriesindx)-effb0vals(1);
        % cf = exp(-effb0val*D);

        % % Normalize and append to Y array
        % thisb0img = p(1)*b0img + p(2);
        % IMGS(seriesindx,:,:,:) = (bimg./thisb0img);
        % % IMGS(seriesindx,:,:,:) = (bimg./thisb0img)./cf;
        
        % IMGS(seriesindx,:,:,:) = (1/p(1))*((1./b0img)-p(2)); 


    end

    clear DINFOS
    clear ImageArrays


    % == Make sample mask
    
    switch samplename
    
        case '20250224_UQ4'
    
            % Cylinder centred at (128, 114)
            
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-128).^2 + (Ys-114).^2 <100^2, 1, 1, szbase(3));
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;
    
        case '20250407_UQ5'
    
            % Cylinder centred at (123, 119), radius 75
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-122).^2 + (Ys-119).^2 <100^2, 1, 1, szbase(3));
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-40:end)=false;
    
    
    
        case '20250414_UQ6'
    
            % Cylinder centred at (122, 117), radius 75
    
            samplemask = zeros(szbase);
            
            [Xs, Ys] = meshgrid( ...
               1:szbase(2), ...
               1:szbase(1) ...
                );
    
            samplemask(:,:,:) = repmat((Xs-122).^2 + (Ys-117).^2 <100^2, 1, 1, szbase(3));
            samplemask(:,:,1:20)=false;
            samplemask(:,:,end-10:end)=false;
    end




    % == Construct composition map
    
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
    
    thiscomposition = reshape(COMPOSITION,[prod(szmap), 3]);
    bool =  sum(thiscomposition,2)>0;
    thiscomposition = thiscomposition(bool, :);
    
    thisimgs = reshape(IMGS, [Nimg, prod(szimg)]);
    thisimgs = thisimgs(:,bool);
       
    composition = [composition; thiscomposition];
    imgs = [imgs, thisimgs];

end


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
Dl = 2.1e-3;
err = 0.001;
% Initialise array for signal measurements
signals = zeros(3, Nimg, 3);

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
    % [beta_fit, resnorm, residual, exitflag, output, lambda, jacobian]=lsqcurvefit(func, beta0, X, y, lb, ub);
    % beta_fit = mdl.Coefficients.Estimate;


    y=y-X(:,3)*Sl;
    mdl = fitlm(X(:,1:2), y, 'Intercept', false);
    beta_fit = mdl.Coefficients.Estimate;
    signals(1:2,imgindx,1) = beta_fit;
    signals(3,imgindx,1) = Sl;


    % % Test predicted signals
    % y_pred = X(:,3)*Sl + sum( X(:,1:2).*repmat(beta_fit', size(X,1),1) ,2);
    % f=figure;
    % scatter(y,y_pred)
    % ylim([0 1])
    % xlim([0 1])
    % corr(y,y_pred)
    % close(f);

    % % Standard error calculation
    % mse = resnorm / (length(y) - length(beta0));
    % cov_matrix = mse * inv(jacobian' * jacobian);
    % stderr = sqrt(diag(cov_matrix));  
    % 
    % signals(:,imgindx,2) = stderr;

    
    % BOOTSTRAPPING TEST
    N=length(y);
    B=5001;
    bootstrap_indices = randi(N, N, B);
    BootFits = zeros(B,3);

    for bindx = 1:B
        % [thisbeta_fit]=lsqcurvefit(func, beta0, X(bootstrap_indices(:,bindx), :), y(bootstrap_indices(:,bindx)), lb, ub);
        % BootFits(bindx, :) = thisbeta_fit;

        thismdl = fitlm(X(bootstrap_indices(:,bindx), 1:2), y(bootstrap_indices(:,bindx)), 'Intercept', false);
        thisbeta_fit = thismdl.Coefficients.Estimate;
        BootFits(bindx, 1:2) = thisbeta_fit;
        BootFits(bindx, 3) = Sl;
    end


    signals(:,imgindx,2) = std(BootFits);
    signals(:,imgindx,3) = prctile(BootFits,95)-prctile(BootFits,5);

    %
    % % TESTING NEURAL NETWORK
    % net = fitnet(10);
    % net = train(net, X', y');
    % y_pred = net(X');
    % figure;
    % scatter(y, y_pred');
    % 
    % 
    % % test E signal
    % EX = zeros(size(X));
    % EX(:,2)=linspace(0,1,size(X,1));
    % 
    % y_pred = net(X');
    % 
    % figure
    % scatter(linspace(0,1,size(X,1)), y_pred);




end



%% Display measured signals

f=figure;
bshift = 5;

% Short Delta

indices = 2:6;
s = signals(:,indices,:);
bvals = [scheme(indices).bval];

errorbar(bvals-1*bshift, s(1,:,1), s(1,:,3), '-*', color=[0.4660 0.6740 0.1880], DisplayName='S (Short \Delta)');
hold on
errorbar(bvals-1*bshift, s(2,:,1), s(2,:,3), '-*', color=[0.8500 0.3250 0.0980], DisplayName='G (Short \Delta)');
errorbar(bvals-1*bshift, s(3,:,1), s(3,:,3), '-*', color=[0 0.4470 0.7410], DisplayName = 'L (Short \Delta)');

% Long Delta
indices = 7:11;
s = signals(:,indices,:);

errorbar(bvals+1*bshift, s(1,:,1), s(1,:,3), '-.*', color=[0.4660 0.6740 0.1880], DisplayName='S (Long \Delta)');
errorbar(bvals+1*bshift, s(2,:,1), s(2,:,3), '-.*', color=[0.8500 0.3250 0.0980], DisplayName='G (Long \Delta)');
errorbar(bvals+1*bshift, s(3,:,1), s(3,:,3), '-.*', color=[0 0.4470 0.7410], DisplayName = 'L (Long \Delta)');

xticks(bvals); 
xticklabels(bvals)
ylim([-0.05,0.7])
ylabel('Normalized signal estimate')
xlim([800,2200])
xlabel('b-value (s/mm^{2})')
xticks(bvals)
xticklabels(["1000", "1250", "1500", "1750", "2000"])
grid on
legend('NumColumns', 2);
ax=gca();
ax.FontSize=12;
f.Position = [488   242   720   480];


%% Save measured signals

savefolder = fullfile(projectfolder, "Outputs", "Signal Measurement", schemename);
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