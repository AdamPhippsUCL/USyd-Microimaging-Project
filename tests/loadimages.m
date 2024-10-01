% MATLAB script to load images

% DICOM folder
dfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\USyd Microimaging Project\Imaging Data\20240711_141030_RB_MouseProstate_d20240711_1_1\3\pdata\1\dicom";

% All DICOM filenames in folder
x = dir(dfolder);
fnames = cell2mat(transpose({char(x(3:end).name)}));
dfolders = repmat([char(dfolder) '/' ], [length(fnames), 1]);
dfnames = [string([dfolders  fnames])];
dfnames = cellstr(transpose(dfnames));

dinfo = dfparse(dfnames);

img = d2mat(dinfo, {'slice','TE'}, 'op','dv');
TEvec = [dinfo([dinfo.sl] == 1).EchoTime];

figure
imshow(img(:,:,25,1), [])