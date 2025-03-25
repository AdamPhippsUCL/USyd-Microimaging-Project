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
    % % Long DELTA
    % 'SE_b0_SPOIL10%',...
    % 'STEAM_LongDELTA_40',...
    % 'STEAM_LongDELTA_60',...
    % 'STEAM_LongDELTA_80',...
    % 'STEAM_LongDELTA_100',...
    % 'STEAM_LongDELTA_120'
    ...
    % % Short DELTA
    % 'SE_b0_SPOIL10% (640 micron)',...
    % 'STEAM_ShortDELTA_15 (640 micron)',...
    % 'STEAM_ShortDELTA_20 (640 micron)',...
    % 'STEAM_ShortDELTA_30 (640 micron)',...
    % 'STEAM_ShortDELTA_40 (640 micron)',...
    % 'STEAM_ShortDELTA_50 (640 micron)'
    ...

    % All DELTA
    'SE_b0_SPOIL10% (640 micron)',...
    'STEAM_ShortDELTA_15 (640 micron)',...
    'STEAM_ShortDELTA_20 (640 micron)',...
    'STEAM_ShortDELTA_30 (640 micron)',...
    'STEAM_ShortDELTA_40 (640 micron)',...
    'STEAM_ShortDELTA_50 (640 micron)',...
    'STEAM_LongDELTA_40 (640 micron)',...
    'STEAM_LongDELTA_60 (640 micron)',...
    'STEAM_LongDELTA_80 (640 micron)',...
    'STEAM_LongDELTA_100 (640 micron)',...
    'STEAM_LongDELTA_120 (640 micron)'...
    ...
    % % Mixed DELTA 1
    % 'SE_b0_SPOIL10% (1600 micron)',...
    % 'STEAM_ShortDELTA_15 (1600 micron)',...
    % 'STEAM_ShortDELTA_30 (1600 micron)',...
    % 'STEAM_ShortDELTA_50 (1600 micron)',...
    % 'STEAM_LongDELTA_60 (1600 micron)',...
    % 'STEAM_LongDELTA_100 (1600 micron)',...
        ...
    % % Mixed DELTA 2
    % 'SE_b0_SPOIL10% (1600 micron)',...
    % 'STEAM_ShortDELTA_20 (1600 micron)',...
    % 'STEAM_ShortDELTA_40 (1600 micron)',...
    % 'STEAM_LongDELTA_40 (1600 micron)',...
    % 'STEAM_LongDELTA_80 (1600 micron)',...
    % 'STEAM_LongDELTA_120 (1600 micron)'...

    %  % Mixed DELTA 3
    % 'SE_b0_SPOIL10% (1600 micron)',...
    % 'STEAM_ShortDELTA_40 (1600 micron)',...
    % 'STEAM_ShortDELTA_50 (1600 micron)',...
    % 'STEAM_LongDELTA_40 (1600 micron)',...
    % 'STEAM_LongDELTA_60 (1600 micron)',...
    % 'STEAM_LongDELTA_80 (1600 micron)'...
    };

% Denoised data
UseDenoisedData = true;


%% Processing details

% Model type
% modeltype = 'ADC';
modeltype = 'RDI - 1 compartment - 2 param';


% Scheme name
schemename = '20250224_UQ4 AllDELTA';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";
load(fullfile(schemesfolder, schemename));

% Fitting technique
fittingtechnique = 'MLP';

% Model folder
modelsfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models');
% switch modeltype
%     case {'RDI - 1 compartment - 2 param','RDI - 2 compartment - 4 param', 'ADC'}
%         modelsfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models');
%     case {'MFP v2'}
%         modelsfolder = fullfile(projectfolder, 'Scripts', 'MFP', fittingtechnique, 'models');
% end



% Average diffusion directions?
averagedirections = true;
Ndirec = 6;


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
switch averagedirections

    case true

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
        % b0img = imgaussfilt(b0img, 0.5);



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


    case false

        % Initialise Y matrix
        Y = ones([Ndirec size(ImageArrays(1).ImageArray, 1:3), length(scheme)]);
        bvec = zeros(1, 2*(length(SeriesDescriptions)-1));
        
        % seriesindx = 1 used for data normalisation!
        img = ImageArrays(1).ImageArray;
        dinfo = DINFOS(1).dinfo;
        bvals = [dinfo(:).DiffusionBValue];
        b0bools = (bvals==0);
        b0imgs = img(:,:,:,b0bools);
        b0img = mean(b0imgs,4);
        meanb0 = mean(b0img(:));
        % 
        % % Apply Gaussian smoothing
        % b0img = imgaussfilt(b0img, 0.5);



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
            bval = bvals(find(bbools,1));
            bvec(2*(seriesindx-1))=bval;
            bimgs = img(:,:,:,bbools);


            for direcindx = 1:Ndirec
                
                bimg = bimgs(:,:,:,direcindx);
                % Normalize and append to Y array
                Y(direcindx,:,:,:,2*(seriesindx-1)) = (meanb0/thismeanb0)*(bimg./b0img); 

            end
            
        end
        
        % Check scheme agreement
        if ~all(bvec == [scheme(:).bval])
            error('Scheme does not match data')
        end
        


end


%% Model fitting

% Define MLP model folder
modelfolder = fullfile(modelsfolder, modeltype, schemename);

switch modeltype

    case 'RDI - 1 compartment - 2 param'

        [R, dIC] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );

    case 'RDI - 2 compartment - 4 param'

        [fIC, R, dIC, dEES] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );


    case 'RDI - 2 compartment - 3 param'

        [fIC, R, d] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );
        

    case 'RDI - 3 compartment - fixed params'

        [fs, fg, fl] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );


    case 'ADC'

        switch averagedirections
            
            case true
                [ADC, S0] = calcADC(Y, [scheme(:).bval]);

            case false
                
                ADC = zeros([Ndirec, size(b0img)]);
                for direcindx = 1:Ndirec
                    ADC(direcindx,:,:,:) = calcADC(squeeze(Y(direcindx,:,:,:,:)), [scheme(:).bval]);
                end
    
        end


    case 'No VASC VERDICT (AMICO)'

        [fIC, fEES, fVASC, R] = verdict_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder...
            );


    case {'MFP v1', 'MFP v2'}

        [R, D] = mfp_fit( ...
            Y, ...
            scheme,...
            modeltype = modeltype,...
            fittingtechnique=fittingtechnique,...
            modelfolder = modelfolder...
            );



end


%% Save results

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting');

outf = fullfile(outputfolder, SampleName, modeltype, schemename, fittingtechnique);
mkdir(outf)


switch modeltype

    case 'RDI - 1 compartment - 2 param'

        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');        

    case 'RDI - 2 compartment - 4 param'
        
        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');
        save(fullfile(outf, 'dEES.mat'), 'dEES');

        % Cellularity
        C = fIC./(R.^3);
        save(fullfile(outf, 'C.mat'), 'C');


    case 'RDI - 2 compartment - 3 param'
        
        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'd.mat'), 'd');

        % Cellularity
        C = fIC./(R.^3);
        save(fullfile(outf, 'C.mat'), 'C');   

    case 'RDI - 3 compartment - fixed params'

        save(fullfile(outf, 'fs.mat'), 'fs');
        save(fullfile(outf, 'fg.mat'), 'fg');
        save(fullfile(outf, 'fl.mat'), 'fl');

    case 'ADC'

        save(fullfile(outf, 'ADC.mat'), 'ADC');
        save(fullfile(outf, 'S0.mat'), 'S0');

    case 'No VASC VERDICT (AMICO)'
        
        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'fEES.mat'), 'fEES');
        save(fullfile(outf, 'fVASC.mat'), 'fVASC');
        save(fullfile(outf, 'R.mat'), 'R');


    case {'MFP v1', 'MFP v2'}

        save(fullfile(outf, 'D.mat'), 'D');
        save(fullfile(outf, 'R.mat'), 'R');

end