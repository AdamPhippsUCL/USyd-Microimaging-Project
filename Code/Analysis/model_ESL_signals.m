% Script to apply modelling to estimated ESL signals

clear;
projectfolder = pwd;

%% Initialise results structure

RESULTS = struct();

SaveRESULTS = true;

%% Sample and scheme details

outputfolder = fullfile(projectfolder, 'Outputs', 'ESL Signal Estimation'); 

samplename = 'Multi-sample'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6';

signals = load(fullfile(outputfolder, samplename, 'signals.mat')).signals;
scheme = load(fullfile(outputfolder, samplename, 'scheme.mat')).scheme;
nscheme = length(scheme);

%% Modelling details

components = {'E', 'S'};

modelnames = {
  'ADC',...
  'DKI',...
  'Sphere',...
  'Ball+Sphere'...
    };

% Fitting
lambda = 0e-3; % Regularisation
fittingtechnique = 'LSQ';

DisplayPredictions = false;

for compindx = 1:length(components)

    component = components{compindx};

    switch component
        case 'E'
            indx=1;
        case 'S'
            indx=2;
        case 'L'
            indx=3;
    end
    
    s = signals(indx, :, 1);
    s_LCI = signals(indx, :, 3);
    s_UCI = signals(indx, :, 4);
    err = signals(indx, :, 2);

    if DisplayPredictions

        f=figure;
        bshift=25;
        switch component
            case 'E'
                color = [0.8500 0.3250 0.0980];
                T = 'Epithelium';
            case 'S'
                color = [0.4660 0.6740 0.1880];
                T = 'Stroma';
            case 'L'
                color = [0    0.4470    0.7410];
                T = 'Lumen';
        end
        
        errorbar( ...
            [scheme(2:6).bval]-0*bshift, ...
            ...log(s(2:6)), ...
            ...log(s(2:6))-log(s_LCI(2:6)), ...
            ...log(s_UCI(2:6))-log(s(2:6)), ...
            (s(2:6)), ...
            (s(2:6))-(s_LCI(2:6)), ...
            (s_UCI(2:6))-(s(2:6)), ...
            '-*', ...
            color=color, MarkerSize=6, ...
            DisplayName = 'Estimated (Short \Delta)')

        hold on

        errorbar( ...
            [scheme(7:end).bval]-0*bshift, ...
            ...log(s(7:end)),...
            ...log(s(7:end))-log(s_LCI(7:end)), ...
            ...log(s_UCI(7:end))-log(s(7:end)), ...
            (s(7:end)),...
            (s(7:end))-(s_LCI(7:end)), ...
            (s_UCI(7:end))-(s(7:end)), ...
            '--*', color=color, MarkerSize=6, ...
            DisplayName = 'Estimated (Long \Delta)')
        
    end

    for modindx = 1:length(modelnames)

        modelname = modelnames{modindx};

        % == Initial guess and bounds
        switch modelname     
            case 'ADC'            
                Nparam = 2;
                beta0 = [1, 1];
                lb = [0,0];
                ub=[2,3];       
                
            case 'DKI'        
                Nparam = 3;
                beta0 = [1, 1, 0.5];
                lb = [0, 0, 0];
                ub = [2, 3, 5];

            case 'Sphere'
                Nparam = 3;
                beta0 = [20, 0.5, 1];
                lb = [1, 0, 0.8];
                ub = [100, 3, 1.2];

            case 'Ball+Sphere'
                Nparam = 5;
                beta0 = [0.1, 10, 1, 1, 1];
                lb = [0, 1, 0, 0, 0.8];
                ub = [1, 50, 3, 3, 1.2];
        end

        
        % == Model fitting 
        
        % Modelling predictions
        [params, resnorm] = fitting_func( ...
            s, ...
            scheme, ...
            modelname = modelname, ...
            fittingtechnique = fittingtechnique,...
            Nparam=Nparam,...
            beta0=beta0,...
            lb=lb,...
            ub=ub,...
            lambda=lambda ...
            );
        
        % AIC
        AIC = nscheme*log(resnorm/nscheme) + 2*Nparam;

        % == Error estimation 
        
        func = @(x) fitting_func( ...
            x,...
            scheme, ...
            modelname=modelname,...
            fittingtechnique=fittingtechnique,...
            lb=lb,...
            ub=ub,...
            Nparam=Nparam,...
            beta0=beta0,...
            lambda=lambda);
        
        % Define input step and bounds
        step = 0.01*ones(size(s));
        xlb = zeros(size(s));
        xub = ones(size(s));
        
        % Estimate Jacobian at signals
        J = JacobianEst(func, s, step=step, xlb=xlb, xub=xub);
        
        % Estimate parameter errors
        signal_var = diag(err.^2);
        params_var =  J*(signal_var*J');
        params_err = sqrt(diag(params_var));

        % == Format results
        
        n = length(RESULTS)+1;
        if ~numel(fieldnames(RESULTS))
            n = 1;
        end
        RESULTS(n).SampleName = samplename;
        RESULTS(n).Component = component;
        RESULTS(n).ModelName = modelname;
        RESULTS(n).ModelParams = params;
        RESULTS(n).ParamError = transpose(params_err);
        RESULTS(n).FitResidual = resnorm;
        RESULTS(n).AIC=AIC;

        if DisplayPredictions

            % Predicted vs measured signals
            pred = zeros(size(signals(indx, :, 1)));

            for ischeme = 2:nscheme
                
                bval = scheme(ischeme).bval;
                delta = scheme(ischeme).delta;
                DELTA = scheme(ischeme).DELTA;
    
                pred(ischeme) = diffusion_model(params, [bval, delta, DELTA], modelname = modelname);
    
            end
    
            switch modelname
                case 'ADC' 
                    markercolor =  [0.5, 0.5, 0.5]	;
                    marker = 'o';
                    markersize = 20;
                    thisbshift = 1*bshift;
                case 'DKI'
                    ...
                case 'Sphere'
                    ...
                case 'Ball+Sphere'
                    markercolor = 'k';%[0.3010, 0.7450, 0.9330];
                    marker = 'x';
                    markersize = 50;
                    thisbshift = 1*bshift;
            end

            scatter( ...
                [scheme(2:end).bval]+thisbshift, ...
                ...log(pred(2:end)), 50,'x', ...
                (pred(2:end)), markersize,marker, ...
                MarkerEdgeColor=markercolor, ...
                LineWidth=1.5, ...
                DisplayName = ['Predicted (' modelname ')'])
            title(T)
            legend
            xticks([1000, 1250, 1500, 1750, 2000])
            xlim([900, 2100])
            xlabel('b-value (s/mm^2)')
            ylabel('dMRI signal')
            switch component
                case 'E'
                    ylim([0.39, 0.66])
                case 'S'
                    ylim([0.22, 0.52])
            end
            grid on
            ax = gca();
            ax.FontSize = 12;
            saveas(f, fullfile(projectfolder, 'Figures', ['Model_Predictions_' T '_.png']))

        end

        % OLD CODE NOT UPDATED

        % % == Profile likelihood
        % 
        % fixedindx = 1;
        % 
        % Nvals = 100;
        % low = (40*lb(fixedindx)+ub(fixedindx))/40;
        % high = (40*ub(fixedindx)+lb(fixedindx))/40;
        % fixedvals = linspace(low,high,Nvals);%linspace(lb(fixedindx), ub(fixedindx),Nvals);
        % valspacing = fixedvals(2)-fixedvals(1);
        % resnorms = zeros(1,Nvals);
        % 
        % for findx = 1:Nvals
        % 
        %     fixedval = fixedvals(findx);
        % 
        %     thisbeta0 = beta0;
        %     thisbeta0(fixedindx) = fixedval;
        % 
        %     thislb = lb;
        %     thislb(fixedindx) = fixedval-0.5*valspacing;
        % 
        %     thisub = ub;
        %     thisub(fixedindx) = fixedval+0.5*valspacing;
        % 
        %     % Modelling predictions
        %     [params, resnorm] = fitting_func( ...
        %         s, ...
        %         scheme, ...
        %         modeltype = modeltype, ...
        %         fittingtechnique = fittingtechnique,...
        %         beta0=thisbeta0,...
        %         lb=thislb,...
        %         ub=thisub,...
        %         lambda=lambda ...
        %         );
        % 
        %     resnorms(findx)=resnorm;
        % 
        % end
        %
        % figure
        % scatter(fixedvals, resnorms)
        % ylim([0, 0.005])
        % %   

    end
end

% Save RESULTS
if SaveRESULTS
    folder = fullfile(projectfolder, 'Outputs', 'Model Fitting', 'ESL signal profiles', samplename);
    mkdir(folder);
    save(fullfile(folder, 'RESULTS.mat'), 'RESULTS')
end


%% Model fitting function (for Jacobian error estimation)

function [params, resnorm] = fitting_func(signals, scheme, opts)

    arguments
        signals 
        scheme
        opts.modelname
        opts.fittingtechnique
        opts.Nparam
        opts.beta0
        opts.lb
        opts.ub
        opts.lambda = 0
    end

    Nparam = length(opts.beta0);

    Y = reshape(signals, [1,1,1,length(signals)]);

    [outputs{1:Nparam}, resnorm, ] = diffusion_model_fit( ...
        Y, ...
        scheme, ...
        modelname=opts.modelname, ...
        fittingtechnique=opts.fittingtechnique,...
        Nparam = opts.Nparam,...
        beta0=opts.beta0,...
        lb=opts.lb,...
        ub=opts.ub,...
        lambda=opts.lambda);

    params = cell2mat(outputs);

end



