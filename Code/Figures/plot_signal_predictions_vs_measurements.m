% Display error on signal predictions (using ESL) compared to signal
% measurements (For benign samples)

clear;
projectfolder = pwd;

%% Sample and image details

% Samples
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'};

% Sample groups
Cancer = {'4B', '4M',  '6N' };
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };

% Image
seriesindx =7;
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
CancerPred = [];
CancerMeasure = [];
CancerCOMP = [];

BenignPred = [];
BenignMeasure = [];
BenignCOMP = [];

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

        if ismember(SampleShortNum, Cancer)
            disp(SampleShortNum)
            CancerPred = [CancerPred; this_pred];
            CancerMeasure = [CancerMeasure; this_measure];
            CancerCOMP = [CancerCOMP; this_COMP];
        
        elseif ismember(SampleShortNum, Benign) 
            BenignPred = [BenignPred; this_pred];
            BenignMeasure = [BenignMeasure; this_measure];
            BenignCOMP = [BenignCOMP; this_COMP];

        end

    end

end


%% MAKE FIGURES

% Benign

f=figure;

scatter(BenignPred, BenignMeasure,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=  BenignCOMP);
hold on
plot([0 0.8],[0, 0.8], color = [.1 .1 .1], LineStyle = '--', LineWidth = 1.2);

ylim([-0.02, 0.8])
xlim([-0.02, 0.8])
grid on
xlabel('Predicted Signal')
ylabel('Measured Signal')
title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])

text(0.025, 0.97, ['R^2 = ' sprintf( '%0.3f', R2(1)) ' (' sprintf('%0.3f', R2(2)) ', ' sprintf('%0.3f', R2(3)) ')'], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

f.Position = [488   242   660   400];
ax = gca();
ax.FontSize = 12;

saveas(f, fullfile(projectfolder, 'Figures', ['Predicted vs Measured signal b' num2str(bval) '_Delta' num2str(DELTA) '.png']))


%% Bland Altman

avg = (BenignPred+BenignMeasure)/2;
diff = (BenignMeasure-BenignPred);

% Limits of agreement
upperLOA = mean(diff)+1.96*std(diff);
lowerLOA = mean(diff)-1.96*std(diff);

f=figure;
scatter(avg, diff ,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=  BenignCOMP, HandleVisibility='off');
yline(mean(diff), '-', DisplayName='Bias', LineWidth=1)
hold on
yline(mean(diff)+1.96*std(diff), '-.', DisplayName='95% LOA', LineWidth=1)
yline(mean(diff)-1.96*std(diff), '-.', HandleVisibility="off", LineWidth=1)
xlim([0, 0.8])
ylim([-0.4, 0.4])
legend
grid on
xlabel('Mean of Predicted and Measured Signal')
ylabel('Measured Signal - Predicted Signal ')
title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
ax = gca();
ax.FontSize = 12;
f.Position = [488   242   660   400];
saveas(f, fullfile(projectfolder, 'Figures', ['Signal Bland Altman b' num2str(bval) '_Delta' num2str(DELTA) '.png']))
