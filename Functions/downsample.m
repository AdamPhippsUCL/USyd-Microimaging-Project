% MATLAB function for downsampling an image

function imgout = downsample(img, window, opts)


arguments
        
    img
    window
    opts.overlap = true;
    opts.filter = 'None'
    opts.sigma = 1;

end

% img = [nx, ny, nz, nscheme]
% window = [nx, ny, nz] 

% Make sure window defined in 3D
if length(window) ~= 3
    window = [window 1];
end

% If multiple image series inputted
if length(size(img))==4
    nscheme = size(img, 4);
else
    nscheme = 0;
end


% Get image size
[Nx, Ny, Nz] = size(img, 1:3);

% Get window size
wx = window(1);
wy = window(2);
wz = window(3);


% Define filter
switch opts.filter

    case 'None'

        filter = ones(wx, wy, wz);


    case 'Gaussian'

        [x,y,z] = ndgrid(-(wx-1)/2:(wx-1)/2, -(wy-1)/2:(wy-1)/2, -(wz-1)/2:(wz-1)/2);
        d = x.^2 + y.^2 + z.^2;
        filter = exp(-d/(2*opts.sigma^2));



end

sumfilt = sum(filter(:));


% Loop and downsample

switch opts.overlap

    case false

        % Calculate number of window positions
        nx = floor(Nx/wx);
        ny = floor(Ny/wy);
        nz = floor(Nz/wz);
        
        if nscheme == 0
            imgout = zeros(nx, ny, nz);
        else
            imgout = zeros(nx, ny, nz, nscheme);
        end


        for indx = 1:nx
            for jndx = 1:ny
                indx
                jndx
                for kndx = 1:nz
                    
                    if nscheme == 0
                        vals = filter.*img((indx-1)*wx+1:(indx)*wx, (jndx-1)*wy+1:(jndx)*wy, (kndx-1)*wz+1:(kndx)*wz);
                        val = sum(vals(:))/sumfilt;
                        imgout(indx, jndx, kndx) = val;
                    else
                        for schndx = 1:nscheme
                            vals = filter.*img((indx-1)*wx+1:(indx)*wx, (jndx-1)*wy+1:(jndx)*wy, (kndx-1)*wz+1:(kndx)*wz, schndx);
                            val = sum(vals(:))/sumfilt;
                            imgout(indx, jndx, kndx, schndx) = val;
                        end
                    end
        
        
                end
            end
        end


    case true

        if nscheme == 0
            imgout = zeros(Nx-wx+1, Ny-wy+1, Nz-wz+1);
        else
            imgout = zeros(Nx-wx+1, Ny-wy+1, Nz-wz+1, nscheme);
        end


        for indx = 1:Nx-wx+1
            for jndx = 1:Ny-wy+1
                for kndx = 1:Nz-wz+1

                    if nscheme == 0
                        vals = filter.*img(indx:indx+wx-1, jndx:jndx+wy-1, kndx:kndx+wz-1);
                        val = sum(vals(:))/sumfilt;
                        imgout(indx, jndx, kndx) = val;       
                    else
                        for schndx = 1:nschemes
                            vals = filter.*img(indx:indx+wx-1, jndx:jndx+wy-1, kndx:kndx+wz-1, schndx);
                            val = sum(vals(:))/sumfilt;
                            imgout(indx, jndx, kndx, schndx) = val;    
                        end
                    end

                    


                end
            end
        end






end