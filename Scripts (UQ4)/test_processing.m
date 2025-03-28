% TEST SCRIPT to implement idea discussed with Roger

clear;

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

samplename = '20250224_UQ4';

%% Base image

% (Image upon which ROIs are defined)
seriesdescription = '3DMGE_20u';

% Load image
baseImg = load(fullfile(projectfolder, 'Imaging Data', 'MAT', samplename, seriesdescription, 'avgImageArray.mat')).avgImageArray;

szbase = size(baseImg);


%% Load masks

maskfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data\MAT\20250224_UQ4\3DMGE_20u\masks\SAMPLE1";

samplemask = load(fullfile(maskfolder, 'sample.mat')).img;
glandularmask = load(fullfile(maskfolder, 'glandular.mat')).img;
mediummask = load(fullfile(maskfolder, 'medium.mat')).img;

% Medium signal value
mediumsignal = prctile(baseImg(logical(mediummask)), 1);

% Define updated masks
GLANDULAR = and( logical(glandularmask), baseImg<mediumsignal);
GLANDULAR = and(GLANDULAR, logical(samplemask));
LUMEN = and(baseImg>mediumsignal, logical(samplemask));
STROMA = and( ~logical(GLANDULAR), ~logical(LUMEN));

% displaymasks = zeros([size(GLANDULAR), 3]);
% displaymasks(:,:,:,1) = logical(GLANDULAR);
% displaymasks(:,:,:,2) = logical(STROMA);
% displaymasks(:,:,:,3) = logical(LUMEN);

%% Load parameter maps

% Model type
modeltype = 'ADC';

% Scheme name
schemename = '20250224_UQ4 LongDELTA';

% Fitting technique
fittingtechnique = 'LSQ';

% Parameter folder
paramfolder = fullfile(projectfolder, "Outputs", "Model Fitting", samplename, modeltype, schemename, fittingtechnique);

% Load parameter map
parammap = load(fullfile(paramfolder, 'ADC.mat')).ADC;


szmap = size(parammap);

%% Find low-res voxels composition

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



%% TEST SCATTER

sl = 15:17;

% Parameters
paramvals = parammap(sl,:,:);
% fICvals = fIC(sl,:,:);
% Rvals = R(sl,:,:);
% dICvals = dIC(sl,:,:);
% dEESvals = dEES(sl,:,:);

mask = squeeze(sum(COMPOSITION(sl,:,:,:),4)>0);
STROMAvals = COMPOSITION(sl,:,:,1);
GLANDULARvals = COMPOSITION(sl,:,:,2);
LUMENvals = COMPOSITION(sl,:,:,3);

figure
scatter([paramvals(mask)], [GLANDULARvals(mask)], '*', HandleVisibility = 'off')
hold on
% xlabel('ADC (s/mm^2)')
xlabel('fIC')
ylabel('Glandular tissue fraction')
ylim([min([GLANDULARvals(mask)]) - 0.05*range([GLANDULARvals(mask)]), max([GLANDULARvals(mask)]) + 0.05*range([GLANDULARvals(mask)])])


% figure
% scatter([fICvals(mask)], [GLANDULARvals(mask)], '*', DisplayName = 'Glandular')
% hold on
% % scatter([fICvals(mask)], [STROMAvals(mask)], '*', DisplayName = 'Stroma')
% % hold on
% % scatter([fICvals(mask)], [LUMENvals(mask)], '*', DisplayName = 'Lumen')
% xlabel('fIC')
% legend;
% 
% 
% 
% figure
% scatter([Rvals(mask)], [GLANDULARvals(mask)], '*', DisplayName = 'Glandular')
% hold on
% % scatter([Rvals(mask)], [STROMAvals(mask)], '*', DisplayName = 'Stroma')
% % hold on
% % scatter([Rvals(mask)], [LUMENvals(mask)], '*', DisplayName = 'Lumen')
% xlabel('R')
% legend;

% 
% figure
% scatter([dICvals(mask)], [GLANDULARvals(mask)], '*', DisplayName = 'Glandular')
% hold on
% scatter([dICvals(mask)], [STROMAvals(mask)], '*', DisplayName = 'Stroma')
% hold on
% scatter([dICvals(mask)], [LUMENvals(mask)], '*', DisplayName = 'Lumen')
% xlabel('dIC')
% legend;
% 
% 
% figure
% scatter([dEESvals(mask)], [GLANDULARvals(mask)], '*', DisplayName = 'Glandular')
% hold on
% scatter([dEESvals(mask)], [STROMAvals(mask)], '*', DisplayName = 'Stroma')
% hold on
% scatter([dEESvals(mask)], [LUMENvals(mask)], '*', DisplayName = 'Lumen')
% xlabel('dEES')
% legend;

%% Correleation and LOBF

% LOBF
p = polyfit([paramvals(mask)], [GLANDULARvals(mask)], 1);
vals = linspace( min([paramvals(mask)]) - 0.05*range([paramvals(mask)]), max([paramvals(mask)]) + 0.05*range([paramvals(mask)]) , 2);
plot(vals, p(1)*vals+p(2), DisplayName='Linear Regression');

% Correlation coefficient
[R,P, RL, RU] = corrcoef([paramvals(mask)], [GLANDULARvals(mask)]);
R = R(1,2);
P = P(1,2);
RL = RL(1,2);
RU = RU(1,2);

legendHandle = legend(Location='northeast');


% Get legend position
legendPos = legendHandle.Position; 

% Define textbox position directly below legend
textboxPos = [legendPos(1)-0.01, legendPos(2) - legendPos(4) + 0.00, legendPos(3) + 0.1, legendPos(4)];

annotation('textbox', textboxPos, 'String', ['CC = ' sprintf('%.3f', R) ' (' sprintf('%.3f', RL) ', ' sprintf('%.3f', RU) ')'], ...
    'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontSize', 10);