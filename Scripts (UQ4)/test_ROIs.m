% TEST script to extract model parameter estimate values from ROIs

clear;

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

samplename = '20250224_UQ4';

%% Base image

% (Image upon which ROIs are defined)
seriesdescription = '3DMGE_20u (40 micron)';

% Load image
baseImg = load(fullfile(projectfolder, 'Imaging Data', 'MAT', samplename, seriesdescription, 'avgImageArray.mat')).avgImageArray;

szbase = size(baseImg);


%% Load parameter maps

% Model type
modeltype = 'RDI - 2 compartment - 4 param';

% Scheme name
schemename = '20250224_UQ4 LongDELTA';

% Fitting technique
fittingtechnique = 'MLP';

% Parameter folder
paramfolder = fullfile(projectfolder, "Outputs", "Model Fitting", samplename, modeltype, schemename, fittingtechnique);

% Load parameter maps
fIC = load(fullfile(paramfolder, 'fIC.mat')).fIC;
R = load(fullfile(paramfolder, 'R.mat')).R;
dIC = load(fullfile(paramfolder, 'dIC.mat')).dIC;
dEES = load(fullfile(paramfolder, 'dEES.mat')).dEES;

szmap = size(fIC);


%% ROI

% ROI name
ROIname = 'STROMA1';

% Load ROI
ROI = load(fullfile(projectfolder, 'Imaging Data', 'MAT', samplename, seriesdescription, [ROIname '.mat'])).img;

sum(ROI(:))

%% Find low-res voxels contained by ROI.


ROImask = zeros(szmap);
ResFactor = szbase./szmap;

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
        for slindx = 1:szmap(3)

            baserows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
            basecols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
            baseslices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));


            ROIvals = ROI(baserows, basecols, baseslices);

            if all(ROIvals(:)==1)
                ROImask(rindx, cindx, slindx)=1;
            end
            


        end
    end
end







