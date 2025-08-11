% Function to apply RDI fitting

function varargout = diffusion_model_fit(Y, scheme, opts)

% Adam Phipps, rmapajp@ucl.ac.uk

arguments

    Y % Normalised signal data array [nx, ny, nz, nscheme]
    scheme % scheme structure

    % ===== OPTIONS

    % Fitting
    opts.modelname
    opts.fittingtechnique  
    opts.mask = []

    % LSQ fitting
    opts.Nparam
    opts.beta0
    opts.lb
    opts.ub
    opts.lambda % regularisation
    opts.calcstderr = false;



end

% Remove NaNs and Infs from Y
Y(isnan(Y)) = 0;
Y(isinf(Y)) = 0;

% Read scheme, data array, and image sizes
nscheme = length(scheme) ;
szY = size(Y) ;
szmap = szY(1:end-1) ;

if ~exist('opts','var') || isempty(opts.mask)
    opts.mask = ones(szmap) ;
end

try
    Y = Y.*opts.mask;
catch
    disp('Mask size mismatch')
end

% Reshape data array (flatten each image)
Y = reshape(Y,[prod(szmap) nscheme]) ;

%% Fitting

switch opts.fittingtechnique

    case 'LSQ'

        % Define predictors
        X = zeros(nscheme, 3);
        X(:,1) = [scheme(:).bval];
        X(:,2) = [scheme(:).delta];
        X(:,3) = [scheme(:).DELTA];

        options = optimset('Display', 'off');

        Nparam = opts.Nparam;

        func = @(b,X) diffusion_model_func(b, X, opts.modelname, @diffusion_model, opts.beta0, opts.lambda);

        % Initialise array for results
        N = prod(szmap);
        x = zeros(N, Nparam);
        stderr = zeros(N, Nparam);
        RESNORM = zeros(N, 1);

        for indx = 1:N

            % disp(indx/N)
            
            thisy = transpose([Y(indx, :), zeros(1, Nparam)]);
            [thisx, resnorm, residual, ~, ~, ~, jacobian] = lsqcurvefit(func, opts.beta0, X, thisy, opts.lb, opts.ub, options);
            x(indx,:)=thisx;
            RESNORM(indx)=resnorm;

            if opts.calcstderr
                dof = length(residual) - length(thisx);
                mse = resnorm / dof;
                cov_matrix = mse * inv(jacobian' * jacobian);
                stderr(indx,:) = sqrt(diag(cov_matrix));
            end

        end  

        % Read results
        for n = 1:Nparam
            varargout{n} = reshape(x(:,n),szmap);
        end
        varargout{Nparam+1} = reshape(RESNORM(:,1), szmap);
        
        if opts.calcstderr
            for n = 1:Nparam
                varargout{Nparam+1+n} = reshape(stderr(:,n), szmap);
            end
        end

end


end


% Define model function for least squares optimization
function Y_out = diffusion_model_func(b,X, modelname, diffusion_model, beta0, lambda)  

    N = size(X,1);
    Y_out = zeros([N,1]);

    for i=1:N
        Y_out(i) = diffusion_model(b,X(i,:),modelname=modelname);
    end

    % Regularisation
    Y_out = [Y_out; transpose(lambda.*((b-beta0)./beta0))];
end
