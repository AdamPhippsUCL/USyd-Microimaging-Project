% Figure for showing correlation coefficients for linear model

clear;
projectfolder = pwd;

SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6'};

% Output folder
outputfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement');

for sindx = 1:length(SampleNames)
    
    samplename = SampleNames{sindx};

    % Load signal measurement results
    thisRESULTS = load(fullfile(outputfolder, samplename, 'RESULTS.mat' )).RESULTS;

    thisR2s=arrayfun(@(s) str2num(s.R2(1:6)), thisRESULTS(2:end)); 

    if sindx == 1
        R2s = thisR2s;
    else
        R2s = [R2s; thisR2s];
    end


end

figure
boxplot(R2s' )
ylim([0, 1])
ylabel('R^2')
xticklabels(SampleNames)