% Script to measure S, G, and L signals

clear;
projectfolder = pwd;

%% Sample and image details

% Sample
multisample = true;
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'};
% multisample = false;
% SampleNames = {'20250524_UQ9'};

schemename = '20250224_UQ4 AllDELTA';
schemesfolder = fullfile(projectfolder, 'Schemes');
load(fullfile(schemesfolder, schemename));
nscheme = length(scheme);
Nimg = nscheme;

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
samplenums = [];

for sindx = 1:length(SampleNames)

    samplename = SampleNames{sindx};

    % == Load masks

    baseseriesdescription = '3DMGE_20u';
    baseimg = load(fullfile(ImagingDataFolder, samplename, baseseriesdescription, 'avgImageArray.mat' )).avgImageArray;

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

        normimg = (1/p(1))*(bimg./(b0img));
        IMGS(seriesindx,:,:,:) =  normimg;

        % Save normalized image
        ImageArray = normimg;
        save(fullfile(ImagingDataFolder, samplename, SeriesDescriptions{seriesindx}, 'normalisedImageArray.mat'), 'ImageArray')


    end

    clear DINFOS
    clear ImageArrays


    % == Make sample mask

    [Ys, Xs, Zs] = ndgrid( ...
        1:szbase(2),...
        1:szbase(1),...
        1:szbase(3)...
        );

    switch samplename
    
        case '20250224_UQ4'
    
            % Cylinder centred at (128, 114), radius 75 (1.5mm)

            samplemask = (Xs-128).^2 + (Ys-114).^2 <75^2 ;
        
            % Remove ends
            % samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;

            % === Define masks for individual samples

            Bmask = samplemask.*and(Zs>=1, Zs<208);
            Mmask = samplemask.*and(Zs>=208, Zs<398);
            Nmask = samplemask.*and(Zs>=398, Zs<640);

    
        case '20250407_UQ5'
    
            % Cylinder centred at (123, 119), radius 75 (1.5mm)

            samplemask = (Xs-123).^2 + (Ys-119).^2 <75^2 ; 

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,600:end)=false;


            % == Define masks for individual samples

            Bmask = samplemask.*and(Zs>=1, Zs<136);
            Mmask = samplemask.*and(Zs>=136, Zs<344);
            Nmask = samplemask.*and(Zs>=344, Zs<641);
    
    
    
        case '20250414_UQ6'
    
            % Cylinder centred at (122, 117), radius 76 (1.52mm)

            samplemask = (Xs-122).^2 + (Ys-117).^2 <76^2;

            % Remove ends
            samplemask(:,:,1:20)=false;
            % samplemask(:,:,end-10:end)=false;

            % == Define masks for individual samples

            Bmask = samplemask.*and(Zs>=15, Zs<285);
            Mmask = samplemask.*and(Zs>=285, Zs<552);
            Nmask = samplemask.*and(Zs>=552, Zs<641);
    
    
        case '20250522_UQ7'
    
            % Cylinder centred at (124, 117)

            samplemask = (Xs-124).^2 + (Ys-117).^2 <75^2;

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;

            % EXCLUDE REGION DUE TO BUBBLES
            rows = 1:90;
            cols = 46:198;
            slices = 225:510;
            samplemask(rows, cols, slices)=false;

            % == Define masks for individual samples

            Bmask = samplemask.*and(Zs>=1, Zs<228);
            Mmask = samplemask.*and(Zs>=228, Zs<506);
            Nmask = samplemask.*and(Zs>=506, Zs<641);            


        case '20250523_UQ8'
    
            % Cylinder centred at (123, 119)

            samplemask = (Xs-123).^2 + (Ys-119).^2 <75^2;

            % Remove ends
            samplemask(:,:,1:10)=false;
            samplemask(:,:,end-10:end)=false;            


            % == Define masks for individual samples

            Bmask = samplemask.*and(Zs>=1, Zs<198);
            Mmask = samplemask.*and(Zs>=198, Zs<458);
            Nmask = samplemask.*and(Zs>=458, Zs<641);   

        case '20250524_UQ9'
    
            % Cylinder centred at (118, 127), radius 75 (1.5mm)
    
            samplemask = (Xs-118).^2 + (Ys-127).^2 <80^2;

            % REMOVE TOP AND BOTTOM REGIONS OF MEDIUM
            samplemask(:,:,1:40)=false;
            samplemask(:,:,570:end)=false;   

            % == Define masks for individual samples

            Bmask = samplemask.*and(Zs>=48, Zs<330);
            Mmask = samplemask.*and(Zs>=330, Zs<330);
            Nmask = samplemask.*and(Zs>=330, Zs<562);   

    end

    % == Construct composition map
    
    % Map size
    szmap = szimg;
    COMPOSITION = zeros([szmap, 3]);
    ResFactor = szbase./szmap;
    Nvoxel = prod(ResFactor);

    NMASK = zeros(szmap);
    MMASK = zeros(szmap);
    BMASK = zeros(szmap);

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

                % MASKs for B, M, N samples
                thisNmask = Nmask(baserows, basecols, baseslices);
                thisMmask = Mmask(baserows, basecols, baseslices);
                thisBmask = Bmask(baserows, basecols, baseslices);

                if sum(thisNmask(:)) > 0.5*Nvoxel
                    NMASK(rindx, cindx, slindx) = 1;
                elseif sum(thisMmask(:)) > 0.5*Nvoxel
                    MMASK(rindx, cindx, slindx) = 1;
                elseif sum(thisBmask(:)) > 0.5*Nvoxel
                    BMASK(rindx, cindx, slindx) = 1;
                end
    
            end
        end
    end


    % Save COMPOSITION
    folder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, SeriesDescriptions{1});
    mkdir(folder)
    save(fullfile(folder, 'COMPOSITION.mat'), 'COMPOSITION');

    % Save sample masks
    save(fullfile(folder, 'NMASK.mat'), 'NMASK');
    save(fullfile(folder, 'MMASK.mat'), 'MMASK');
    save(fullfile(folder, 'BMASK.mat'), 'BMASK');




    % NEW
    % Select voxels IN BENIGN SAMPLES AND WITH NON-ZERO composition

    % Cancer samples = UQ4B, UQ4M, UQ6N
    switch samplename(end-2:end)
        case 'UQ4'
            % Only include sample N
            COMPOSITION = COMPOSITION.*double(NMASK);
        case 'UQ6'
            % Only include samples B and M
            COMPOSITION = COMPOSITION.*(double(BMASK)+double(MMASK));
    end


    % Select voxels with non-zero composition
    thiscomposition = reshape(COMPOSITION,[prod(szmap), 3]);
    bool =  sum(thiscomposition,2)>0;
    thiscomposition = thiscomposition(bool, :);

    thisimgs = reshape(IMGS, [Nimg, prod(szimg)]);
    thisimgs = thisimgs(:,bool);


    % Sample masks
    this_bmask = BMASK(:);
    this_bmask = this_bmask(bool);
    this_mmask = MMASK(:);
    this_mmask = this_mmask(bool);
    this_nmask = NMASK(:);
    this_nmask = this_nmask(bool);

    samplenums = [samplenums; (3*sindx-2)*this_bmask + (3*sindx-1)*this_mmask + (3*sindx)*this_nmask];
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
Dl = 2.0e-3;
err = 0.001;

% Initialise array for signal measurements
signals = zeros(3, Nimg, 4);

% Initialise RESULTS structure
RESULTS = struct();

for imgindx = 1:Nimg

    RESULTS(imgindx).SampleName = samplename;
    RESULTS(imgindx).bval = scheme(imgindx).bval;
    RESULTS(imgindx).DELTA = scheme(imgindx).DELTA;

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

    % Test linear model assumption
    R2 = mdl.Rsquared.Ordinary;
    residuals = mdl.Residuals.Raw;

    % BLAND ALTMAN
    measured = y+X(:,3)*Sl;
    predicted = y-residuals+X(:,3)*Sl;
    avg = (measured+predicted)/2;
    diff = (measured-predicted);
    LOA = [mean(diff), mean(diff)-1.96*std(diff), mean(diff)+1.96*std(diff)];

    f=figure;
    scatter(predicted, measured)
    hold on
    plot(y,y)
    close(f);

    signals(1:2,imgindx,1) = beta_fit;
    signals(3,imgindx,1) = Sl;
    
    % BOOTSTRAPPING 

    % All data
    N=length(y);

    % % Samples
    % N = max(samplenums);

    B=100;
    bootstrap_indices = randi(N, N, B);
    BootFits = zeros(B,3);
    BootR2s = zeros(B,1);


    for bindx = 1:B

        % [thisbeta_fit]=lsqcurvefit(func, beta0, X(bootstrap_indices(:,bindx), :), y(bootstrap_indices(:,bindx)), lb, ub);
        % BootFits(bindx, :) = thisbeta_fit;

        thismdl = fitlm(X(bootstrap_indices(:,bindx), 1:2), y(bootstrap_indices(:,bindx)), 'Intercept', false);

        % this_indices = bootstrap_indices(:,bindx);
        % ids = [];
        % for i = 1:length(this_indices)
        %     indx = this_indices(i);
        %     ids = [ids; find(samplenums == indx)];
        % 
        % end
        % thismdl = fitlm(X(ids, 1:2), y(ids), 'Intercept', false);


        thisbeta_fit = thismdl.Coefficients.Estimate;
        BootFits(bindx, 1:2) = thisbeta_fit;
        BootFits(bindx, 3) = Sl;

        % thisy_pred = X(:,3)*Sl + sum( X(:,1:2).*repmat(thisbeta_fit', size(X,1),1) ,2);
        BootR2s(bindx)= thismdl.Rsquared.Ordinary;

    end


    signals(:,imgindx,2) = std(BootFits); % Standard error
    signals(:,imgindx,3) = prctile(BootFits,1); % 2.5th percentile
    signals(:,imgindx,4) = prctile(BootFits,99); % 97.5th percentile

    % RESULTS
    RESULTS(imgindx).E = [signals(2, imgindx, 1), signals(2, imgindx, 3), signals(2, imgindx, 4)];% [num2str(signals(2, imgindx, 1)) ' (' num2str(signals(2, imgindx, 2)) ')'];
    RESULTS(imgindx).S = [signals(1, imgindx, 1), signals(1, imgindx, 3), signals(1, imgindx, 4)];% [num2str(signals(1, imgindx, 1)) ' (' num2str(signals(1, imgindx, 2)) ')'];
    RESULTS(imgindx).L = [signals(3, imgindx, 1), signals(3, imgindx, 3), signals(3, imgindx, 4)];% [num2str(signals(3, imgindx, 1)) ' (' num2str(signals(3, imgindx, 2)) ')'];
    RESULTS(imgindx).R2 = [R2, prctile(BootR2s, 2.5), prctile(BootR2s, 97.5)];%[num2str(R2) ' (' num2str(std(BootR2s)) ')'];
    RESULTS(imgindx).Residuals = residuals;
    RESULTS(imgindx).y = y;

    % LOA
    RESULTS(imgindx).LOA = LOA;



end

RESULTS(1).X = X;


%% Display measured signals

f=figure;
bshift=5;

% Short Delta

indices = 2:6;
s = signals(:,indices,:);
bvals = [scheme(indices).bval];

% Display error bars with 95% confidence intervals from botstrapping
errorbar(bvals-1*bshift, s(1,:,1), s(1,:,1)-s(1,:,3), s(1,:,4)-s(1,:,1), '-*', color=[0.4660 0.6740 0.1880], DisplayName='S (Short \Delta)');
hold on
errorbar(bvals-1*bshift, s(2,:,1), s(2,:,1)-s(2,:,3), s(2,:,4)-s(2,:,1), '-*', color=[0.8500 0.3250 0.0980], DisplayName='G (Short \Delta)');
errorbar(bvals-1*bshift, s(3,:,1), s(3,:,1)-s(3,:,3), s(3,:,4)-s(3,:,1),'-*', color=[0 0.4470 0.7410], DisplayName = 'L (Short \Delta)');


% % Display error bars with STANDARD ERROR from botstrapping
% errorbar(bvals-1*bshift, s(1,:,1), s(1,:,2), '-*', color=[0.4660 0.6740 0.1880], DisplayName='S (Short \Delta)');
% hold on
% errorbar(bvals-1*bshift, s(2,:,1), s(2,:,2), '-*', color=[0.8500 0.3250 0.0980], DisplayName='G (Short \Delta)');
% errorbar(bvals-1*bshift, s(3,:,1), s(3,:,2),'-*', color=[0 0.4470 0.7410], DisplayName = 'L (Short \Delta)');


% Long Delta
indices = 7:11;
s = signals(:,indices,:);

% Display error bars with 95% confidence intervals from botstrapping
errorbar(bvals+1*bshift, s(1,:,1), s(1,:,1)-s(1,:,3), s(1,:,4)-s(1,:,1), '-.*', color=[0.4660 0.6740 0.1880], DisplayName='S (Long \Delta)');
errorbar(bvals+1*bshift, s(2,:,1), s(2,:,1)-s(2,:,3), s(2,:,4)-s(2,:,1),  '-.*', color=[0.8500 0.3250 0.0980], DisplayName='G (Long \Delta)');
errorbar(bvals+1*bshift, s(3,:,1), s(3,:,1)-s(3,:,3), s(3,:,4)-s(3,:,1), '-.*', color=[0 0.4470 0.7410], DisplayName = 'L (Long \Delta)');

% % Display error bars with STANDARD ERROR from botstrapping
% errorbar(bvals+1*bshift, s(1,:,1), s(1,:,2), '-.*', color=[0.4660 0.6740 0.1880], DisplayName='S (Long \Delta)');
% errorbar(bvals+1*bshift, s(2,:,1), s(2,:,2),  '-.*', color=[0.8500 0.3250 0.0980], DisplayName='G (Long \Delta)');
% errorbar(bvals+1*bshift, s(3,:,1), s(3,:,2), '-.*', color=[0 0.4470 0.7410], DisplayName = 'L (Long \Delta)');

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

switch multisample
    case false
        samplename = SampleNames{1};
        savefolder = fullfile(projectfolder, "Outputs", "Signal Measurement", samplename);
        mkdir(savefolder);
    case true
        savefolder = fullfile(projectfolder, "Outputs", "Signal Measurement", 'Multi-sample');
        mkdir(savefolder);
end

% Meta info
Meta = struct();
Meta.SampleNames = SampleNames;
Meta.SeriesDescriptions = SeriesDescriptions;
Meta.Nboot = B;

save(fullfile(savefolder, 'signals.mat'), 'signals')
save(fullfile(savefolder, 'scheme.mat'), 'scheme')
save(fullfile(savefolder, 'Meta.mat'), 'Meta')
save(fullfile(savefolder, 'RESULTS.mat'), 'RESULTS')


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