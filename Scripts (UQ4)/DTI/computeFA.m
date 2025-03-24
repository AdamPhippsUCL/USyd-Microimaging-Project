function FA = computeFA(D_tensor)
    % Compute Fractional Anisotropy (FA) for a 3D volume of diffusion tensors
    %
    % Input:
    %   - D_tensor: 4D array [Nx, Ny, Nz, 6] containing diffusion tensor elements
    %               (Dxx, Dyy, Dzz, Dxy, Dxz, Dyz) for each voxel.
    %
    % Output:
    %   - FA: 3D array [Nx, Ny, Nz] of FA values.

    [Nx, Ny, Nz, ~] = size(D_tensor);
    FA = zeros(Nx, Ny, Nz);  % Preallocate FA map

    % Loop over each voxel
    for x = 1:Nx
        for y = 1:Ny
            for z = 1:Nz
                D_voxel = squeeze(D_tensor(x, y, z, :));
                
                if any(D_voxel)  % Check for valid tensor values
                    % Reconstruct 3x3 diffusion tensor
                    D = [D_voxel(1), D_voxel(4), D_voxel(5);
                         D_voxel(4), D_voxel(2), D_voxel(6);
                         D_voxel(5), D_voxel(6), D_voxel(3)];
                    
                    % Compute eigenvalues of the tensor
                    lambda = eig(D);
                    
                    % Ensure positive eigenvalues (to avoid numerical issues)
                    if all(lambda > 0)
                        MD = mean(lambda);  % Mean diffusivity
                        FA(x, y, z) = sqrt(3/2) * sqrt(sum((lambda - MD).^2)) / sqrt(sum(lambda.^2));
                    end
                end
            end
        end
    end
end