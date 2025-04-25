% Load D and downsample to size of downsampled images

clear;
projectfolder = pwd;

% Load D
SampleName = '20250224_UQ4';

load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', '40u_DtiSE_2012_SPOIL10% (20 micron)', 'D.mat'));

downsamplewindow = [3*8,3*8,4*8]; 
overlap = false;

D_DS = imresize3(D, [10,10,20]);

% Save
folder = fullfile(projectfolder, 'Scripts', 'SGL Signals', 'D', 'SampleName');
mkdir(folder);
save(fullfile(folder, 'D_DS.mat' ), 'D_DS');