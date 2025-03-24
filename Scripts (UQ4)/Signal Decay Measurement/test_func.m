% FUNCTION for non-linear model fitting

function signals = test_func(b, X)


% Tissue modelled with three compartments
%
% 1. Glandular -> restricted diffusion, radius Rg, diffusivity Dg
% 2. Stroma -> restricted diffusion, radius Rs, diffusivity Ds
% 3. Lumen -> Gaussian diffusion, diffusivity Dl

Rs = b(1);
Ds = b(2);
Rg = b(3);
Dg = b(4);
Rl = b(5);
Dl = b(6);

Nobvs = size(X,1);

signals = zeros(Nobvs,1);

for indx = 1:Nobvs

    bval = X(indx, 1);
    delta = X(indx,2);
    DELTA = X(indx,3);
    G = stejskal(delta, DELTA, bval=bval);

    fs = X(indx,4);
    fg = X(indx,5);
    fl = X(indx,6);

    % STROMA SIGNAL
    Ss = sphereGPD(delta, DELTA, G, Rs, Ds);


    % GLANDULAR SIGNAL
    Sg = sphereGPD(delta, DELTA, G, Rg, Dg);

    % LUMEN SIGNAL
    Sl = sphereGPD(delta, DELTA, G, Rl, Dl);

    % Total signal
    S = fs*Ss + fg*Sg + fl*Sl;

    signals(indx) = S;


end


end