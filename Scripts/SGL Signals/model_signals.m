% Script to apply modelling to measured signals

clear;
projectfolder = pwd;

%% Initialise results structure

RESULTS = struct();


%% Sample and scheme details

outputfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement'); 

samplename = 'Multi-sample'; % 'ESMRMB', '20250224_UQ4', '20250407_UQ5', '20250414_UQ6';


signals = load(fullfile(outputfolder, samplename, 'signals.mat')).signals;
scheme = load(fullfile(outputfolder, samplename, 'scheme.mat')).scheme;
nscheme = length(scheme);

%% Modelling details

components = {'G', 'S'};

modeltypes = {
    'ADC',...
    'DKI',...
    'RDI - 1 compartment - 2 param (S0)',...
    'RDI - 2 compartment - 3 param (S0)',...
    'RDI - 2 compartment - 4 param (S0)'
    };


% modeltype = 'ADC';
% modeltype = 'DKI';
% modeltype = 'RDI - 1 compartment - 2 param (S0)';
% modeltype = 'RDI - 2 compartment - 3 param (S0)';
% modeltype = 'RDI - 2 compartment - 4 param (S0)';

% Fitting
lambda = 0e-3;
fittingtechnique = 'LSQ';


for compindx = 1:length(components)

    component = components{compindx};

    for modindx = 1:length(modeltypes)

        modeltype = modeltypes{modindx};

        % == Initial guess and bounds
        switch modeltype     
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
            case 'RDI - 1 compartment - 2 param'        
                Nparam = 2;
                beta0 = [20, 0.5];
                lb = [1, 0];
                ub = [100, 3];        
            case 'RDI - 2 compartment - 3 param'        
                Nparam = 2;
                beta0 = [0.1, 10, 1];
                lb = [0, 1, 0];
                ub = [1, 50, 3];        
            case 'RDI - 2 compartment - 4 param'       
                Nparam = 4;
                beta0 = [0.1, 10, 1, 1];
                lb = [0, 1, 0, 0];
                ub = [1, 50, 3, 3];               
            case 'RDI - 1 compartment - 2 param (S0)'
                Nparam = 2;
                beta0 = [20, 0.5, 1];
                lb = [1, 0, 0.8];
                ub = [100, 3, 1.2];
            case 'RDI - 2 compartment - 3 param (S0)'
                Nparam = 2;
                beta0 = [0.1, 10, 1,1];
                lb = [0, 1, 0,0.8];
                ub = [1, 50, 3,1.2];
            case 'RDI - 2 compartment - 4 param (S0)'
                Nparam = 5;
                beta0 = [0.1, 10, 1, 1, 1];
                lb = [0, 1, 0, 0, 0.8];
                ub = [1, 50, 3, 3, 1.2];
        end

        
        % == Model fitting 
        
        switch component
            case 'S'
                indx=1;
            case 'G'
                indx=2;
            case 'L'
                indx=3;
        end
        
        s = signals(indx, :, 1);
        err = signals(indx, :, 2);
        
        % Modelling predictions
        [params, resnorm] = fitting_func( ...
            s, ...
            scheme, ...
            modeltype = modeltype, ...
            fittingtechnique = fittingtechnique,...
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
            modeltype=modeltype,...
            fittingtechnique=fittingtechnique,...
            lb=lb,...
            ub=ub,...
            beta0=beta0,...
            lambda=lambda);
        
        % Define input step and bounds
        step = 0.01*ones(size(s));
        xlb = zeros(size(s));
        xub = ones(size(s));
        
        % Estimate Jacobian at signals
        J = JacobianEst(func, s, step=step, xlb=xlb, xub=xub);
        
        % % Covariance matrix
        % CoV = inv(J'*J) * resnorm / (size(s,2) - Nparam);
        
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
        RESULTS(n).ModelType = modeltype;
        RESULTS(n).ModelParams = params;
        RESULTS(n).ParamError = transpose(params_err);
        RESULTS(n).FitResidual = resnorm;
        RESULTS(n).AIC=AIC;

% disp(RESULTS);

        % % == Profile likelihood
        % 
        % fixedindx = 3;
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
        % 
        % figure
        % scatter(fixedvals, resnorms)
        % ylim([0, 0.005])
        %   

    end
end

% Save RESULTS
folder = fullfile(projectfolder, 'Outputs', 'Signal Measurement', samplename, 'Modelling');
mkdir(folder);
save(fullfile(folder, 'RESULTS.mat'), 'RESULTS')




%% Model fitting function (for Jacobian error estimation)

function [params, resnorm] = fitting_func(signals, scheme, opts)

    arguments
        signals 
        scheme
        opts.modeltype
        opts.fittingtechnique
        opts.beta0
        opts.lb
        opts.ub
        opts.lambda = 0
    end

    Nparam = length(opts.beta0);

    Y = reshape(signals, [1,1,1,length(signals)]);

    [outputs{1:Nparam}, resnorm, ] = RDI_fit( ...
        Y, ...
        scheme, ...
        modeltype=opts.modeltype, ...
        fittingtechnique=opts.fittingtechnique,...
        beta0=opts.beta0,...
        lb=opts.lb,...
        ub=opts.ub,...
        lambda=opts.lambda);

    params = cell2mat(outputs);

end



