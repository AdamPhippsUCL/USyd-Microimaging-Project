% Display error on signal predictions (using ESL) compared to signal
% measurements

% Separate into cancer + benign groups

clear;
projectfolder = pwd;


%% Sample and image details

% Samples
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'};

% Sample groups
Cancer = {'4B', '4M',  '6N' };
Benign = {'4N', '5B', '5M', '5N', '6B',  '6M', '7M', '7N', '8B', '8M', '8N', '7B', '9B', '9N' };

% Image
seriesindx =11;
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
        signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'signals.mat')).signals;
        signals = squeeze(signals(:,seriesindx,1));

        % R2 value
        RESULTS = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'RESULTS.mat')).RESULTS;
        R2 = RESULTS(seriesindx).R2;

    % case false
    %     signals = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', SampleName, 'signals.mat')).signals;
    %     signals = squeeze(signals(:,seriesindx,1));
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

% Reformat COMP columns for color coding
CancerCOMP(:, [1 2]) = CancerCOMP(:, [2 1]);
BenignCOMP(:, [1 2]) = BenignCOMP(:, [2 1]);

%% MAKE FIGURES
% 
% % CANCER
% 
% figure
% scatter(CancerPred, CancerMeasure, 6, 'filled', 'MarkerFaceAlpha', 0.7, CData=  CancerCOMP);
% hold on
% plot([0, 0.7],[0, 0.7], '--', color = [0.2, 0.2, 0.2]);
% grid on
% xlim([-0.02, 0.78])
% ylim([-0.02, 0.82])
% xlabel('Predicted Signal (from ESL model)')
% ylabel('Measured Signal')
% title(['Cancer (incl. Gleason \geq 3+3); ', 'b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])

% Benign

f=figure;
scatter(BenignPred, BenignMeasure,  6, 'filled', 'MarkerFaceAlpha', 0.7, CData=  BenignCOMP);
hold on
plot([0 0.8],[0, 0.8], color = [.1 .1 .1], LineStyle = '--', LineWidth = 1.2);
ylim([-0.02, 0.8])
xlim([-0.02, 0.8])
grid on
xlabel('Predicted Signal (from ESL model)')
ylabel('Measured Signal')
title(['b = ' num2str(bval) ' s/mm^2 ; Delta = ' num2str(DELTA) ' ms'])
text(0.025, 0.97, ['R^2 = ' sprintf( '%0.3f', R2(1)) ' (' sprintf('%0.3f', R2(2)) ', ' sprintf('%0.3f', R2(3)) ')'], ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');  % Optional border

saveas(f, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', ['LinMod_b' num2str(bval) '_Delta' num2str(DELTA) '.png']))


%% Bland Altman

avg = (BenignPred+BenignMeasure)/2;
diff = (BenignMeasure-BenignPred);

% Limits of agreement
upperLOA = mean(diff)+1.96*std(diff);
lowerLOA = mean(diff)-1.96*std(diff);

figure
yline(mean(diff), '--')
hold on
yline(mean(diff)+1.96*std(diff))
yline(mean(diff)-1.96*std(diff))
xlim([0, 0.8])
ylim([-0.4, 0.4])

cancer_avg = (CancerPred+CancerMeasure)/2;
cancer_diff = (CancerMeasure-CancerPred);
scatter(cancer_avg, cancer_diff ,  '.', MarkerEdgeAlpha=0.4, CData=CancerCOMP)




% 
% % Variable LOA
% windowSize = 40;
% [sorted_avg, idx] = sort(avg);
% sorted_diff = diff(idx);
% 
% Nwindow = length(avg)/windowSize;
% 
% movAVGs = [mean(sorted_avg(1: windowSize/2))];
% movDIFFs = [mean(sorted_diff(1: windowSize/2))];
% movLOAs = [1.96*std(sorted_diff(1: windowSize/2))];
% 
% for windx = 1:Nwindow
%     movAVGs = [movAVGs mean(sorted_avg((windx-1/2)*windowSize+1: windx*windowSize))];
%     movDIFFs = [movDIFFs mean(sorted_diff((windx-1)*windowSize+1: windx*windowSize))];
%     movLOAs = [movLOAs 1.96*std(sorted_diff((windx-1)*windowSize+1: windx*windowSize))];
% end
% 
% 
% movAVGs = [movAVGs mean(sorted_avg(end-windowSize/2: end))];
% movDIFFs = [movDIFFs mean(sorted_diff(end-windowSize/2: end))];
% movLOAs = [movLOAs 1.96*std(sorted_diff(end-windowSize/2: end))];
% 
% % 
% % plot(movAVGs, movDIFFs, 'k')
% % plot(movAVGs, movDIFFs-movLOAs, 'k')
% % plot(movAVGs, movDIFFs+movLOAs, 'k')
% 
% % p = polyfit(movAVGs, movDIFFs,3);
% % movDIFFs = polyval(p, sorted_avg);
% % 
% % p = polyfit(movAVGs, movLOAs,3);
% % movLOAs = polyval(p, sorted_avg);
% % 
% % plot(sorted_avg, movDIFFs, 'k')
% % plot(sorted_avg, movDIFFs-movLOAs, 'k')
% % plot(sorted_avg, movDIFFs+movLOAs, 'k')
% 
% 
% % Savgol smoothing
% movDIFFs = sgolayfilt(movDIFFs,1,25);
% movLOAs = sgolayfilt(movLOAs,1,25);
% 
% movDIFFs = movmean(movDIFFs, 20);
% movLOAs = movmean(movLOAs, 20);
% 
% plot(movAVGs, movDIFFs, 'k')
% plot(movAVGs, movDIFFs-movLOAs, 'k')
% plot(movAVGs, movDIFFs+movLOAs, 'k')
% 
