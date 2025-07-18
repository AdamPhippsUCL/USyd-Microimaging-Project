% Script to format modelling resilts into table

clear;
projectfolder = pwd;

% Load model fitting results
RESULTS = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'Multi-sample', 'Modelling', 'RESULTS.mat')).RESULTS;
N = length(RESULTS);

% Initialise table
T = table( ...
    'Size', ...
    [N, 4], ...
    'VariableTypes', {'string', 'string', 'string', 'string'}, ...
    'VariableNames', {'Component', 'Model Name', 'Estimated Parameters', 'AIC'});

% Fill in table
for indx = 1:N

    component = RESULTS(indx).Component;
    modelname = RESULTS(indx).ModelType;
    params = RESULTS(indx).ModelParams;
    error = RESULTS(indx).ParamError;
    AIC = RESULTS(indx).AIC;

    % Component
    switch component
        case 'G'
            T{indx, 'Component'} = {'Epithelium'};
        case 'S'
            T{indx, 'Component'} = {'Stroma'};
    end


    % Parameters and errors
    switch modelname
        case 'ADC'
            T{indx, 'Model Name'} = {modelname};
            
            T{indx, 'Estimated Parameters'} = {[...
                'S0 = ' sprintf('%.3f', params(1)) ' (' sprintf('%.4f', error(1)) '); ' ...
                'D = ' sprintf('%.3f', params(2)) ' (' sprintf('%.4f', error(2)) ') mm^2/s' ...
            ]};


        case 'DKI'
            T{indx, 'Model Name'} = {modelname};

        T{indx, 'Estimated Parameters'} = {[...
                'S0 = ' sprintf('%.3f', params(1)) ' (' sprintf('%.4f', error(1)) '); ' ...
                'D = ' sprintf('%.3f', params(2)) ' (' sprintf('%.4f', error(2)) ') mm^2/s; ' ...
                'K = ' sprintf('%.3f', params(3)) ' (' sprintf('%.4f', error(3)) ')' ...
            ]};


        case 'RDI - 1 compartment - 2 param (S0)'
            T{indx, 'Model Name'} = {'Sphere'};
            
            T{indx, 'Estimated Parameters'} = {[...
                'S0 = ' sprintf('%.3f', params(3)) ' (' sprintf('%.4f', error(3)) '); ' ...
                'D = ' sprintf('%.3f', params(2)) ' (' sprintf('%.4f', error(2)) ') mm^2/s; ' ...
                'R = ' sprintf('%.3f', params(1)) ' (' sprintf('%.4f', error(1)) ') µm' ...
            ]};


        case 'RDI - 2 compartment - 4 param (S0)'
            T{indx, 'Model Name'} = {'Ball + Sphere'};
    
            T{indx, 'Estimated Parameters'} = {[...
                'S0 = ' sprintf('%.3f', params(5)) ' (' sprintf('%.4f', error(5)) '); ' ...
                'f = ' sprintf('%.3f', params(1)) ' (' sprintf('%.4f', error(1)) '); ' ...
                'D_in = ' sprintf('%.3f', params(3)) ' (' sprintf('%.4f', error(3)) ') mm^2/s; ' ...
                'R = ' sprintf('%.3f', params(2)) ' (' sprintf('%.4f', error(2)) ') µm; ' ...
                'D_out = ' sprintf('%.3f', params(4)) ' (' sprintf('%.4f', error(4)) ') mm^2/s' ...
            ]};

    end

    % AIC
    T{indx, 'AIC'} = {sprintf('%.3f', AIC)};


end


% Save as excel sheet
writetable(T, fullfile(projectfolder, 'Scripts', 'Paper Figures', 'Figures', 'modelling_results.xlsx'));

