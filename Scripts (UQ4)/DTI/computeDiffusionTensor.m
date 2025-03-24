function D_tensor = computeDiffusionTensor(volb, volb0, bvals, bvecs)
    % Estimate diffusion tensor (D) voxel-wise for whole imaging volume
    % Inputs:
    %   - volb: 4D matrix [Nx, Ny, Nz, 6] (6 diffusion-weighted images)
    %   - volb0: 3D matrix [Nx, Ny, Nz] (b=0 image)
    %   - bvals: [6x1] vector of b-values
    %   - bvecs: [6x3] matrix of gradient directions (each row is [gx, gy, gz])
    % Output:
    %   - D_tensor: 4D matrix [Nx, Ny, Nz, 6] (stores Dxx, Dyy, Dzz, Dxy, Dxz, Dyz)

    [Nx, Ny, Nz, Nd] = size(volb);
    
    % Construct B-matrix for least squares
    B = [-bvals .* (bvecs(:,1).^2), ...
         -bvals .* (bvecs(:,2).^2), ...
         -bvals .* (bvecs(:,3).^2), ...
         -2 * bvals .* (bvecs(:,1) .* bvecs(:,2)), ...
         -2 * bvals .* (bvecs(:,1) .* bvecs(:,3)), ...
         -2 * bvals .* (bvecs(:,2) .* bvecs(:,3))];

    % Preallocate output
    D_tensor = zeros(Nx, Ny, Nz, 6);

    % Loop over all voxels
    for x = 1:Nx
        for y = 1:Ny
            for z = 1:Nz
                S0 = volb0(x, y, z);
                if S0 > 0  % Ensure valid signal
                    S = squeeze(volb(x, y, z, :));  % Extract signal for all directions
                    ln_S = log(S / S0);  % Take log signal ratio
                    D_vec = B \ ln_S;  % Solve least squares system
                    D_tensor(x, y, z, :) = D_vec;  % Store tensor elements
                end
            end
        end
    end
end
