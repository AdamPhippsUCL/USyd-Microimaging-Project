% Script to visualise distribution of errors in linear model over sample

clear;
projectfolder = pwd;

%% Sample and image details

% Sample
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6'};


% Image
seriesindx = 11;
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

scheme = load(fullfile(projectfolder, "Schemes", "20250224_UQ4 AllDELTA.mat")).scheme;
bval = scheme(seriesindx).bval;
DELTA = scheme(seriesindx).DELTA;


%% Calculate errors for each sample

multisample = true;

Errors = [];
Samples = [];
Xs = [];

for sindx = 1:length(SampleNames)

    SampleName = SampleNames{sindx};

    ImageFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, SeriesDescription);
    ImageArray = load(fullfile(ImageFolder,'normalisedImageArray.mat')).ImageArray;
    % dinfo = load(fullfile(ImageFolder,'axialdinfo.mat')).dinfo;

    switch multisample
    
        case true
            signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'signals.mat')).signals;
            signals = squeeze(signals(:,seriesindx,1));

        case false

            signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', SampleName, 'signals.mat')).signals;
            signals = squeeze(signals(:,seriesindx,1));
    end
    
    COMPOSITION = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;


    signals = reshape(signals, [1,1,1,3]);
    pred = sum(COMPOSITION.*repmat(signals, [size(COMPOSITION, 1:3)]), 4);
    error = (pred - ImageArray);
    error = error(pred>0);

    Errors = [Errors; error];
    Samples = [Samples; sindx*ones(size(error))];


    % Component fractions
    X = COMPOSITION(repmat(pred, [1,1,1,3])>0);
    X = reshape(X, [length(X)/3, 3]);
    Xnew = X;
    Xnew(:,1)=X(:,2);
    Xnew(:,2)=X(:,1);
    X=Xnew;
    Xs = [Xs; X];


end
            


figure
swarmchart(Samples, Errors, '*',MarkerEdgeAlpha=0.15,  CData=Xs)
hold on
bxplt=boxplot(Errors, Samples);
set(bxplt, 'LineWidth', 1)
ylim([-0.32, 0.32])
yline(0, 'LineStyle', '--')
xlabel('Sample')
ylabel('Predicted signal error')
title(['b-value = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])




