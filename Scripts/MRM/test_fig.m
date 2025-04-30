
% Test script for MRM figure

clear;
projectfolder = pwd;


%% Load modelling results

resultsfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement');
SampleNames = {'20250224_UQ4', '20250407_UQ5', '20250414_UQ6'};

for sindx = 1:length(SampleNames)

    samplename = SampleNames{sindx};

    theseRESULTS = load(fullfile(resultsfolder, samplename, 'Modelling', 'RESULTS.mat' )).RESULTS;

    if sindx == 1
        RESULTS = theseRESULTS;
    else
        RESULTS = [RESULTS, theseRESULTS];
    end
end



%% Plot model results


modeltype = 'RDI - 2 compartment - 3 param (S0)';
Nparam = 3;


for paramindx = 1:Nparam
    
    figure;
    
    for sindx = 1:length(SampleNames)
        
        samplename = SampleNames{sindx};
    
        samplebools = strcmp({RESULTS(:).SampleName}, samplename);
        modelbools = strcmp({RESULTS(:).ModelType}, modeltype);
    
        if sindx == 1
            results = RESULTS(and(samplebools, modelbools));
        else
            results = [results, RESULTS(and(samplebools, modelbools))];
        end
        
    end
    
    
    
    component = 'G';
    componentbools = strcmp({results(:).Component}, component);
    paramvals = arrayfun(@(s) s.ModelParams(paramindx), results(componentbools));
    paramerrs  = arrayfun(@(s) s.ParamError(paramindx), results(componentbools));
    errorbar([1:length(SampleNames)]-0.1, paramvals, paramerrs, 'LineStyle', 'none', 'Marker', '*', DisplayName=component)
    hold on
    
    component = 'S';
    componentbools = strcmp({results(:).Component}, component);
    paramvals = arrayfun(@(s) s.ModelParams(paramindx), results(componentbools));
    paramerrs  = arrayfun(@(s) s.ParamError(paramindx), results(componentbools));
    errorbar([1:length(SampleNames)]+0.1, paramvals, paramerrs, 'LineStyle', 'none', 'Marker', '*',  DisplayName=component)
    hold on
    
    legend
    title(['Param: ' num2str(paramindx)])
    
    xticks([1:length(SampleNames)]);
    xticklabels(SampleNames)

end
