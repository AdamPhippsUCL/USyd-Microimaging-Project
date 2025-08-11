function signal = diffusion_model(model_params, sequence_params, opts)

arguments

   model_params % vector of model parameters (model dependent)
   sequence_params % vector of sequence parameters [bval, delta, DELTA]

   opts.modelname % model type

end

% Unpack sequence parameters
bval = sequence_params(1);
delta = sequence_params(2);
DELTA = sequence_params(3);
G = stejskal(delta, DELTA, bval=bval);

% signal simulation
switch opts.modelname

    case 'ADC'

        % model_params = [S0, D]
        S0 = model_params(1);
        D = model_params(2);

        signal = S0*ball(bval, D);

    case 'DKI'

        % model_params = [S0, D, K]
        S0 = model_params(1);
        D = model_params(2);
        K = model_params(3);

        signal = S0*exp(-bval*D*1e-3 + (1/6)*K*((D*1e-3)^2)*(bval^2));


    case 'Sphere'

        % model_params = [R, D]
        R = model_params(1);
        D = model_params(2);
        S0 = model_params(3);

        signal = S0*sphereGPD(delta, DELTA, G, R, D);


    case 'Ball+Sphere - 3 param'

        % model_params = [fs, R, D]
        fs = model_params(1);
        fb = 1-fs;
        R = model_params(2);
        D = model_params(3);
        S0 = model_params(4);

        % Sphere signal
        Ss = sphereGPD(delta, DELTA, G, R, D);

        % Ball signal
        Sb = ball(bval, D);

        % Total signal
        signal = S0*(fs*Ss + fb*Sb);


    case 'Ball+Sphere'

        % model_params = [S0, fs, R, Ds, Db, S0]
        fs = model_params(1);
        fb = 1-fs;
        R = model_params(2);
        Ds = model_params(3);
        Db = model_params(4);
        S0 = model_params(5);

        % sphere signal
        Ss = sphereGPD(delta, DELTA, G, R, Ds);

        % ball signal
        Sb = ball(bval, Db);

        % Total signal
        signal = S0*(fs*Ss + fb*Sb);


end


end