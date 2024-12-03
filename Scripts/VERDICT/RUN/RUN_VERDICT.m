% MATLAB script to run VERDICT processing

clear;

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20241128_UQ1';

% Series descriptions
SeriesDescriptions = {...
    '40u_verdict_seq1_v2',...
    '40u_verdict_seq2_v2',...
    '40u_verdict_seq3_v2',...
    '40u_verdict_seq4_v2',...
    '40u_verdict_seq5_v2',...
    };


%% VERDICT processing details

verdictfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Scripts\VERDICT";

% Model type
modeltype = 'Original VERDICT';

% Scheme name
schemename = 'UQ Scheme v2';

% Fitting technique
fittingtechnique = 'MLP';

% Model folder
modelsfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Scripts\VERDICT\MLP\models";

%% Data preprocessing

% == Load images and dinfo

ImageArrays = struct();
DINFOS = struct();

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    % Load image and dinfo
    thisfolder = fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription);
    ImageArray = load(fullfile(thisfolder, 'axialImageArray.mat')).ImageArray;
    dinfo = load(fullfile(thisfolder, 'axialdinfo.mat')).dinfo;

    % Append to structures
    DINFOS(seriesindx).dinfo = dinfo;
    ImageArrays(seriesindx).ImageArray = ImageArray;

    
end

% == Load scheme

scheme = load(fullfile(verdictfolder, 'Schemes', [schemename '.mat'])).scheme;


% == Construct Y matrix

% Initialise Y matrix
Y = ones([size(ImageArrays(1).ImageArray, 1:3), length(scheme)]);
bvec = zeros(1, 2*length(SeriesDescriptions));

for seriesindx = 1:length(SeriesDescriptions)

    % Load image
    img = ImageArrays(seriesindx).ImageArray;

    % Load dinfo
    dinfo = DINFOS(seriesindx).dinfo;
    bvals = [dinfo(:).DiffusionBValue];

    % b0 imgs
    b0bools = (bvals==0);
    b0imgs = img(:,:,:,b0bools);
    b0img = mean(b0imgs,4);
    bvec(2*seriesindx-1)=0;

    % b imgs
    bbools = (bvals > 0);
    bimgs = img(:,:,:,bbools);
    bimg = mean(bimgs,4);
    bval = bvals(find(bbools,1));
    bvec(2*seriesindx)=bval;


    % Normalize and append to Y array
    Y(:,:,:,2*seriesindx) = bimg./b0img; 

   
end

% Check scheme agreement
if ~all(bvec == [scheme(:).bval])
    error('Scheme does not match data')
end



%% VERDICT Fitting

% Define MLP model folder
modelfolder = fullfile(modelsfolder, modeltype, schemename);

[fIC, fEES, fVASC, R] = verdict_fit( ...
    Y, ...
    scheme, ...
    modeltype = modeltype,...
    fittingtechnique = fittingtechnique,...
    modelfolder = modelfolder...
    );





