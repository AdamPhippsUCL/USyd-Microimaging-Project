% MATLAB script to run VERDICT processing

clear;

projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250131_CPC1';

% Series descriptions
SeriesDescriptions = {...
    % % STEAM
    % '192700001',...
    % '192720001',...
    % '192670001',...
    % SE
    '192690001',...
    '192710001',...
    '192730001',...
    };

% Denoised data
UseDenoisedData = false;

%% VERDICT processing details

% Model type
modeltype = 'No VASC VERDICT (AMICO)';

% Scheme name
schemename = '20250131_CPC SE';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";
load(fullfile(schemesfolder, schemename));

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
    switch UseDenoisedData
        case true
            thisfolder = fullfile(ImagingDataFolder, 'MAT DN', SampleName, SeriesDescription);
        case false
            thisfolder = fullfile(ImagingDataFolder, 'MAT', SampleName, SeriesDescription);
    end
    ImageArray = load(fullfile(thisfolder, 'axialImageArray.mat')).ImageArray;
    dinfo = load(fullfile(thisfolder, 'axialdinfo.mat')).dinfo;

    % Append to structures
    DINFOS(seriesindx).dinfo = dinfo;
    ImageArrays(seriesindx).ImageArray = ImageArray;

end


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


%% Save outputs

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting');

outf = fullfile(outputfolder, SampleName, modeltype, schemename, fittingtechnique);
mkdir(outf)

save(fullfile(outf, 'fIC.mat'), 'fIC');
save(fullfile(outf, 'fEES.mat'), 'fEES');
save(fullfile(outf, 'fVASC.mat'), 'fVASC');
save(fullfile(outf, 'R.mat'), 'R');

