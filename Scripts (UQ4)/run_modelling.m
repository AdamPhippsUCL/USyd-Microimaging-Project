% MATLAB script to run DW-MRI modelling on UQ4 imaging data

clear;

projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";


%% Image details

% Imaging data folder 
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Sample name
SampleName = '20250224_UQ4';

% Series descriptions
SeriesDescriptions = {...
    % Long DELTA
    % 'SE_b0_SPOIL10% (640 micron)',...
    % 'STEAM_LongDELTA_40 (640 micron)',...
    % 'STEAM_LongDELTA_60 (640 micron)',...
    % 'STEAM_LongDELTA_80 (640 micron)',...
    % 'STEAM_LongDELTA_100 (640 micron)',...
    % 'STEAM_LongDELTA_120 (640 micron)'
    ...
    % % Short DELTA
    % 'SE_b0_SPOIL10% (640 micron)',...
    % 'STEAM_ShortDELTA_15 (640 micron)',...
    % 'STEAM_ShortDELTA_20 (640 micron)',...
    % 'STEAM_ShortDELTA_30 (640 micron)',...
    % 'STEAM_ShortDELTA_40 (640 micron)',...
    % 'STEAM_ShortDELTA_50 (640 micron)'
    ...
    % % All DELTA
    % 'SE_b0_SPOIL10% (640 micron)',...
    % 'STEAM_ShortDELTA_15 (640 micron)',...
    % 'STEAM_ShortDELTA_20 (640 micron)',...
    % 'STEAM_ShortDELTA_30 (640 micron)',...
    % 'STEAM_ShortDELTA_40 (640 micron)',...
    % 'STEAM_ShortDELTA_50 (640 micron)',...
    % 'STEAM_LongDELTA_40 (640 micron)',...
    % 'STEAM_LongDELTA_60 (640 micron)',...
    % 'STEAM_LongDELTA_80 (640 micron)',...
    % 'STEAM_LongDELTA_100 (640 micron)',...
    % 'STEAM_LongDELTA_120 (640 micron)'...
    ...
    % % Mixed DELTA 1
    % 'SE_b0_SPOIL10% ',...
    % 'STEAM_ShortDELTA_15',...
    % 'STEAM_ShortDELTA_30',...
    % 'STEAM_ShortDELTA_50',...
    % 'STEAM_LongDELTA_60',...
    % 'STEAM_LongDELTA_100',...
        ...
    % Mixed DELTA 2
    'SE_b0_SPOIL10%',...
    'STEAM_ShortDELTA_20',...
    'STEAM_ShortDELTA_40',...
    'STEAM_LongDELTA_40',...
    'STEAM_LongDELTA_80',...
    'STEAM_LongDELTA_120'...
    };

% Denoised data
UseDenoisedData = true;


%% Processing details

% Model type
modeltype = 'RDI - 2 compartment - 4 param';

% Scheme name
schemename = '20250224_UQ4 MixedDELTA2';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";
load(fullfile(schemesfolder, schemename));

% Fitting technique
fittingtechnique = 'MLP';

% Model folder
modelsfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models');



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
bvec = zeros(1, 2*(length(SeriesDescriptions)-1));

% seriesindx = 1 used for data normalisation!
img = ImageArrays(1).ImageArray;
dinfo = DINFOS(1).dinfo;
bvals = [dinfo(:).DiffusionBValue];
b0bools = (bvals==0);
b0imgs = img(:,:,:,b0bools);
b0img = mean(b0imgs,4);
meanb0 = mean(b0img(:));

% Apply Gaussian smoothing
b0img = imgaussfilt(b0img, 0.5);

% seriesindx > 1 used for diffusion data
for seriesindx = 2:length(SeriesDescriptions)

    % Load image
    img = ImageArrays(seriesindx).ImageArray;

    % Load dinfo
    dinfo = DINFOS(seriesindx).dinfo;
    bvals = [dinfo(:).DiffusionBValue];

    % b0 imgs
    b0bools = (bvals==0);
    b0imgs = img(:,:,:,b0bools);
    thisb0img = mean(b0imgs,4);
    thismeanb0 = mean(thisb0img(:));

    % b imgs
    bbools = (bvals > 0);
    bimgs = img(:,:,:,bbools);
    bimg = mean(bimgs,4);
    bval = bvals(find(bbools,1));
    bvec(2*(seriesindx-1))=bval;


    % Normalize and append to Y array
    Y(:,:,:,2*(seriesindx-1)) = (meanb0/thismeanb0)*(bimg./b0img); 

end

% Check scheme agreement
if ~all(bvec == [scheme(:).bval])
    error('Scheme does not match data')
end


%% Model fitting

% Define MLP model folder
modelfolder = fullfile(modelsfolder, modeltype, schemename);


switch modeltype

    case {'RDI - 2 compartment - 4 param', 'RDI - 2 compartment - Fixed Ds'}

        [fIC, fEES, R, dIC, dEES] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );



    case 'ADC'

        [ADC, S0] = calcADC(Y, [scheme(:).bval]);


    case 'No VASC VERDICT (AMICO)'

        [fIC, fEES, fVASC, R] = verdict_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );




end


%% Save results

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting');

outf = fullfile(outputfolder, SampleName, modeltype, schemename, fittingtechnique);
mkdir(outf)


switch modeltype

    case {'RDI - 2 compartment - 4 param', 'RDI - 2 compartment - Fixed Ds'}
        
        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'fEES.mat'), 'fEES');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');
        save(fullfile(outf, 'dEES.mat'), 'dEES');

    case 'ADC'

        save(fullfile(outf, 'ADC.mat'), 'ADC');
        save(fullfile(outf, 'S0.mat'), 'S0');

    case 'No VASC VERDICT (AMICO)'
        
        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'fEES.mat'), 'fEES');
        save(fullfile(outf, 'fVASC.mat'), 'fVASC');
        save(fullfile(outf, 'R.mat'), 'R');


end