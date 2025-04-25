% Script to make table 2 (fitting results)

clear;
projectfolder = pwd;

% Load results
RESULTS = load(fullfile(projectfolder, "Outputs", "Signal Measurement", "RESULTS.mat")).RESULTS;
Nresults = length(RESULTS);


%% Make table

varnames = {'Component', 'Model', 'Estimated Parameters', 'AIC'};
vartypes = {'string', 'string', 'string', 'string'};
Ncol = length(varnames);
T = table('Size', [Nresults Ncol], 'VariableTypes', vartypes, 'VariableNames', varnames);

for rindx = 1:Nresults

    component = RESULTS(rindx).Component;
    model = RESULTS(rindx).ModelType;
    params = RESULTS(rindx).ModelParams;
    errs = RESULTS(rindx).ParamError;
    AIC = RESULTS(rindx).AIC;

    T(rindx, "Component") = {component};
    T(rindx, "Model") = {model};
    T(rindx, "AIC") = {sprintf('%.2f', AIC)};

    paramstr = '';
    
    for indx = 1:length(params)
        paramstr = [paramstr sprintf('%.4f' ,params(indx)) ' (' sprintf('%.4f' ,errs(indx)) '); ' ];
    end

    T(rindx, "Estimated Parameters") = {paramstr};
    

end


% Save as excel sheet
writetable(T, fullfile(projectfolder, "Scripts", "Figures", "table.xlsx"));