% Script to calculate predicted signals for each voxel (using ESL segmentations and aggregate ESL signal estimates)

clear;
projectfolder = pwd;

%% Sample and image details

% Samples
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'};

% Image
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
scheme = load(fullfile(projectfolder, "Schemes", "20250224_UQ4 AllDELTA.mat")).scheme;


for seriesindx = 2:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};
    
    bval = scheme(seriesindx).bval;
    DELTA = scheme(seriesindx).DELTA;
    
    % Load signal measurements
    multisample = true;
    switch multisample
        case true
            signals = load(fullfile(projectfolder, 'Outputs', 'ESL signal estimation', 'Multi-sample', 'signals.mat')).signals;
            signals = squeeze(signals(:,seriesindx,1));
    
            % R2 value
            RESULTS = load(fullfile(projectfolder, 'Outputs', 'ESL signal estimation', 'Multi-sample', 'RESULTS.mat')).RESULTS;
            R2 = RESULTS(seriesindx).R2;
    end
    
    % Initialise arrays
    SampleNums = [];
    Predicted = [];
    Measured = [];
    COMP = [];
    
    % Loop over samples
    for sampleindx = 1:length(SampleNames)
    
        SampleName = SampleNames{sampleindx};
        samplenum = SampleName(end);
    
        ImageFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN', SampleName, SeriesDescription);
        ImageArray = load(fullfile(ImageFolder,'normalisedImageArray.mat')).ImageArray;
    
        % Load composition
        COMPOSITION = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', 'COMPOSITION.mat')).COMPOSITION;
        
        % Predicted signals
        signals = reshape(signals, [1,1,1,3]);
        pred = sum(COMPOSITION.*repmat(signals, [size(COMPOSITION, 1:3)]), 4);
    

        % PER SAMPLE
        samplelabels = {'N', 'M', 'B'};
    
        for slabindx = 1:length(samplelabels)
    
            samplelabel = samplelabels{slabindx};
    
            MASK = load(fullfile(projectfolder, 'Outputs', 'Masks', SampleName, 'SE_b0_SPOIL5% (DS)', [samplelabel 'MASK.mat'])).([samplelabel 'MASK']);
    
            % Predicted signal
            this_pred = pred(logical(MASK));
    
            % Measured signal
            this_measure = ImageArray(logical(MASK));
    
            % COMPOSITION
            this_COMP = reshape(COMPOSITION, [], 3);
            flatMASK = MASK(:);
            this_COMP = this_COMP(logical(flatMASK), :);
    
            SampleShortNum = [samplenum samplelabel];
    
            disp(SampleShortNum)
            SampleNums = [SampleNums; repmat({SampleShortNum}, length(this_pred), 1)];
            Predicted = [Predicted; this_pred];
            Measured = [Measured; this_measure];
            COMP = [COMP; this_COMP];
    
        end
    
    end
    
    % Save measured and predicted signals
    folder = fullfile(projectfolder, 'Outputs', 'Signals');
    mkdir(folder)
    save(fullfile(folder, 'SampleNums.mat'), 'SampleNums')
    save(fullfile(folder, 'COMP.mat'), 'COMP')
    
    seriesfolder = fullfile(folder, SeriesDescription);
    mkdir(seriesfolder);
    save(fullfile(seriesfolder, 'Measured.mat'), 'Measured')
    save(fullfile(seriesfolder, 'Predicted.mat'), 'Predicted')

end
