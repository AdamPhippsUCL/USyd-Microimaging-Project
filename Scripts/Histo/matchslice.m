% Script to match histo and MGE slices

clear;
projectfolder = pwd;


%% UQ9N

sample_name = '20250524_UQ9';

% Load image data
ImageArray  = load(fullfile(pwd, 'Imaging Data', 'MAT DN', sample_name, '3DMGE_20u', 'avgImageArray.mat')).avgImageArray;
  
baseseriesdescription = '3DMGE_20u';
maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', sample_name, baseseriesdescription);
GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN = double(load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN);


sl = 546;
slice = ImageArray(:,:,sl);
maxI = max(slice(:));
minI = min(slice(:));

figure
imshow(squeeze(ImageArray(:,:,sl)), [minI, maxI])

% Define slice centre
Cx = 125;
Cy = 114;
Cz = sl;

% Make x, y, z, coordinate matrix  (ImageArray(y, x, z))
[Nx, Ny, Nz] = size(ImageArray);
xs = linspace(-Cx+1, Nx-Cx, Nx);
ys = linspace(-Cy+1, Ny-Cy, Ny);
zs = linspace(-Cz+1, Nz-Cz, Nz);
[Xs, Ys, Zs] = meshgrid(xs, ys, zs);

voxelcoords = cat(4, Xs, Ys, Zs);
voxelcoords = permute(voxelcoords, [4 1 2 3]); % Now size is (3, Nx, Ny, Nz)
voxelcoords = reshape(voxelcoords, 3, []); 


% % T1 = reflect in y
% M = [1, 0, 0; 0, -1, 0; 0, 0, 1];
% invM = inv(M);
% 

% voxelcoords_new = invM * voxelcoords;
% voxelcoords_new = reshape(voxelcoords_new, 3, Nx, Ny, Nz);
% voxelcoords_new = permute(voxelcoords_new, [2 3 4 1]); 
% 
% Xs_new = voxelcoords_new(:,:,:,1);
% Ys_new = voxelcoords_new(:,:,:,2);
% Zs_new = voxelcoords_new(:,:,:,3);
% 
% ImageArray = interp3(Xs, Ys, Zs, ImageArray, Xs_new, Ys_new, Zs_new, "linear");
% 
% figure
% imshow(squeeze(ImageArray (:,:,sl)), [minI, maxI])


% T2 = rotate about z
alpha = -20; % degrees
alpha = alpha*pi/180; % Radians
axis = [0, 0, 1];
M = RotMat(axis, alpha);
invM = inv(M);

voxelcoords_new = invM * voxelcoords;
voxelcoords_new = reshape(voxelcoords_new, 3, Nx, Ny, Nz);
voxelcoords_new = permute(voxelcoords_new, [2 3 4 1]); 

Xs_new = voxelcoords_new(:,:,:,1);
Ys_new = voxelcoords_new(:,:,:,2);
Zs_new = voxelcoords_new(:,:,:,3);

ImageArray = interp3(Xs, Ys, Zs, ImageArray, Xs_new, Ys_new, Zs_new, "linear");
GLANDULAR = interp3(Xs, Ys, Zs, GLANDULAR, Xs_new, Ys_new, Zs_new, "linear");
STROMA = interp3(Xs, Ys, Zs, STROMA, Xs_new, Ys_new, Zs_new, "linear");
LUMEN = interp3(Xs, Ys, Zs, LUMEN, Xs_new, Ys_new, Zs_new, "linear");

figure
imshow(squeeze(ImageArray (:,:,sl)), [minI, maxI])


% T3 rotation in xy plane
alpha = 7.5; % degrees
alpha = alpha*pi/180; % Radians
axis = [sin(30), cos(30), 0];

M = RotMat(axis, alpha);
invM = inv(M);

voxelcoords_new = invM * voxelcoords;
voxelcoords_new = reshape(voxelcoords_new, 3, Nx, Ny, Nz);
voxelcoords_new = permute(voxelcoords_new, [2 3 4 1]); 

Xs_new = voxelcoords_new(:,:,:,1);
Ys_new = voxelcoords_new(:,:,:,2);
Zs_new = voxelcoords_new(:,:,:,3);

ImageArray = interp3(Xs, Ys, Zs, ImageArray, Xs_new, Ys_new, Zs_new, "linear");
GLANDULAR = interp3(Xs, Ys, Zs, GLANDULAR, Xs_new, Ys_new, Zs_new, "linear");
STROMA = interp3(Xs, Ys, Zs, STROMA, Xs_new, Ys_new, Zs_new, "linear");
LUMEN = interp3(Xs, Ys, Zs, LUMEN, Xs_new, Ys_new, Zs_new, "linear");

figure
imshow(squeeze(ImageArray (:,:,sl)), [minI, maxI])




% T4 rotatation
alpha = 2; % degrees
alpha = alpha*pi/180; % Radians
axis = [sin(30), cos(30), 0];

M = RotMat(axis, alpha);
invM = inv(M);

voxelcoords_new = invM * voxelcoords;
voxelcoords_new = reshape(voxelcoords_new, 3, Nx, Ny, Nz);
voxelcoords_new = permute(voxelcoords_new, [2 3 4 1]); 

Xs_new = voxelcoords_new(:,:,:,1);
Ys_new = voxelcoords_new(:,:,:,2);
Zs_new = voxelcoords_new(:,:,:,3);

ImageArray = interp3(Xs, Ys, Zs, ImageArray, Xs_new, Ys_new, Zs_new, "linear");
GLANDULAR = interp3(Xs, Ys, Zs, GLANDULAR, Xs_new, Ys_new, Zs_new, "linear");
STROMA = interp3(Xs, Ys, Zs, STROMA, Xs_new, Ys_new, Zs_new, "linear");
LUMEN = interp3(Xs, Ys, Zs, LUMEN, Xs_new, Ys_new, Zs_new, "linear");


figure
imshow(squeeze(ImageArray (:,:,sl)), [minI, maxI])
hold on

displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = (GLANDULAR);
displaymasks(:,:,:,2) = (STROMA);
displaymasks(:,:,:,3) = (LUMEN);

mask = imshow(squeeze(displaymasks(:,:,sl,:)));
set(mask, 'AlphaData', 0.2)



function R = RotMat(axis, theta)

    % Ensure axis is a column vector
    axis = axis(:);
    axis = axis / norm(axis);  % Normalize the axis vector
    
    x = axis(1);
    y = axis(2);
    z = axis(3);
    
    c = cos(theta);
    s = sin(theta);
    C = 1 - c;
    
    % Rodrigues' rotation formula
    R = [ ...
        c + x^2*C,     x*y*C - z*s, x*z*C + y*s;
        y*x*C + z*s,   c + y^2*C,   y*z*C - x*s;
        z*x*C - y*s,   z*y*C + x*s, c + z^2*C ];
    
end