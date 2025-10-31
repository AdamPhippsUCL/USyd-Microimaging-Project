% Script to run DW-MRI modelling on imaging data (NOT USED IN CURRENT ANALYSIS)

clear;
projectfolder = pwd;

%% Image details

% Sample name
SampleName = '20250524_UQ9'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

% Series descriptions
SeriesDescriptions = {...
    % All DELTA
    'SE_b0_SPOIL5% (DS)',...
    'STEAM_ShortDELTA_15 (DS)',...
    'STEAM_ShortDELTA_20 (DS)',...
    'STEAM_ShortDELTA_30 (DS)',...
    'STEAM_ShortDELTA_40 (DS)',...
    'STEAM_ShortDELTA_50 (DS)',...
    'STEAM_LongDELTA_40 (DS)',...
    'STEAM_LongDELTA_60 (DS)',...
    'STEAM_LongDELTA_80 (DS)',...
    'STEAM_LongDELTA_100 (DS)',...
    'STEAM_LongDELTA_120 (DS)'...
    ...
    };

% Denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end

% Average diffusion directions?
averagedirections = true;
Ndirec = 6;


%% Processing details

% Model type
modelname = 'Ball+Sphere'; % 'ADC', 'DKI', 'Sphere', 'Ball+Sphere'

% Scheme name
schemename = '20250224_UQ4 AllDELTA';
schemesfolder = fullfile(projectfolder, 'Schemes');
load(fullfile(schemesfolder, schemename));
nscheme = length(scheme);

% Fitting technique
fittingtechnique = 'LSQ';

% === LSQ fitting

switch modelname

    case 'ADC'
        
        Nparam=2;

        D=1; Dlb = 0.2; Dub = 2;
        S0=1; S0lb = 0.8; S0ub = 1.2;

        beta0=[S0,D];
        lb=[S0lb,Dlb];
        ub=[S0ub,Dub];

    case 'DKI'
        
        Nparam=3;
        beta0=[S0,D,K];
        lb=[S0lb,dlb,Klb];
        ub=[S0ub,dub,Kub];

    case 'Sphere'

        Nparam = 2;
        beta0 = [R,D];
        lb = [Rlb,Dlb];
        ub = [Rub,Dub];


    case 'Ball+Sphere'

        Nparam = 5;

        % R and Ds free parameters??
        free_R_Ds = false;

        switch free_R_Ds

            case false

                fs = 0.5; fslb = 0; fsub = 1;
                R = 6.5; Rlb = 6.4; Rub = 6.6;
                Ds = 0.55; Dslb = 0.54; Dsub = 0.56;
                Db = 1; Dblb = 0.1; Dbub = 3;
                S0 = 1; S0lb = 0.9; S0ub = 1.1;

            case true

                fs = 0.5; fslb = 0; fsub = 1;
                R = 6.5; Rlb = 1; Rub = 20;
                Ds = 0.01; Dslb = 0.54; Dsub = 3;
                Db = 1; Dblb = 0.1; Dbub = 3;
                S0 = 1; S0lb = 0.9; S0ub = 1.1;

        end

        beta0 = [fs,R,Ds,Db,S0];
        lb = [fslb,Rlb,Dslb,Dblb,S0lb];
        ub = [fsub,Rub,Dsub,Dbub,S0ub];
        
end

% Regularisation
lambda0=0e-2;
lambda = lambda0*ones(1,Nparam);

%% Data preprocessing

% == Load images and dinfo

ImageArrays = struct();
DINFOS = struct();

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    % Load image and dinfo
    thisfolder = fullfile(ImagingDataFolder, SampleName, SeriesDescription);
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
        bvec = zeros(1, length(scheme));
        
        % seriesindx = 1 used for data normalisation!
        img = ImageArrays(1).ImageArray;
        dinfo = DINFOS(1).dinfo;
        bvals = [dinfo(:).DiffusionBValue];
        b0bools = (bvals==0);
        b0imgs = img(:,:,:,b0bools);
        b0img = mean(b0imgs,4);
        % meanb0 = mean(b0img(:));
    
        b0prctiles = [
            prctile(b0img(:), 10),...
            prctile(b0img(:), 20),...
            prctile(b0img(:), 30),...
            prctile(b0img(:), 40),...
            prctile(b0img(:), 50),...
            prctile(b0img(:), 60),...
            prctile(b0img(:), 70),...
            prctile(b0img(:), 80),...
            prctile(b0img(:), 90),...
        ];


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

            thisb0prctiles = [
                prctile(thisb0img(:), 10),...
                prctile(thisb0img(:), 20),...
                prctile(thisb0img(:), 30),...
                prctile(thisb0img(:), 40),...
                prctile(thisb0img(:), 50),...
                prctile(thisb0img(:), 60),...
                prctile(thisb0img(:), 70),...
                prctile(thisb0img(:), 80),...
                prctile(thisb0img(:), 90),...
            ];

            p = polyfit(b0prctiles, thisb0prctiles, 1);

        
            % b imgs
            bbools = (bvals > 0);
            bimgs = img(:,:,:,bbools);
            bimg = mean(bimgs,4);
            bval = bvals(find(bbools,1));
            bvec(seriesindx)=bval;
        
            % Normalize and append to Y array
            Y(:,:,:,seriesindx) = (1/p(1))*(bimg./(b0img));
        
        end
        

    % case false
    % 
    %     % Initialise Y matrix
    %     Y = ones([Ndirec size(ImageArrays(1).ImageArray, 1:3), length(scheme)]);
    %     bvec = zeros(1, length(scheme));
    % 
    %     % seriesindx = 1 used for data normalisation!
    %     img = ImageArrays(1).ImageArray;
    %     dinfo = DINFOS(1).dinfo;
    %     bvals = [dinfo(:).DiffusionBValue];
    %     b0bools = (bvals==0);
    %     b0imgs = img(:,:,:,b0bools);
    %     b0img = mean(b0imgs,4);
    %     meanb0 = mean(b0img(:));
    % 
    %     % % Apply Gaussian smoothing
    %     % b0img = imgaussfilt(b0img, 0.5);
    % 
    % 
    %     % seriesindx > 1 used for diffusion data
    %     for seriesindx = 2:length(SeriesDescriptions)
    % 
    %         % Load image
    %         img = ImageArrays(seriesindx).ImageArray;
    % 
    %         % Load dinfo
    %         dinfo = DINFOS(seriesindx).dinfo;
    %         bvals = [dinfo(:).DiffusionBValue];
    % 
    %         % b0 imgs
    %         b0bools = (bvals==0);
    %         b0imgs = img(:,:,:,b0bools);
    %         thisb0img = mean(b0imgs,4);
    %         thismeanb0 = mean(thisb0img(:));
    % 
    %         % b imgs
    %         bbools = (bvals > 0);
    %         bval = bvals(find(bbools,1));
    %         bvec(seriesindx)=bval;
    %         bimgs = img(:,:,:,bbools);
    % 
    % 
    %         for direcindx = 1:Ndirec
    % 
    %             bimg = bimgs(:,:,:,direcindx);
    %             % Normalize and append to Y array
    %             Y(direcindx,:,:,:,seriesindx) = (meanb0/thismeanb0)*(bimg./b0img); 
    % 
    %         end
    % 
    %     end
    % 
    %     % % Check scheme agreement
    %     % if ~all(bvec == [scheme(:).bval])
    %     %     error('Scheme does not match data')
    %     % end
        
end


%% Model fitting

switch modelname

    case 'ADC'

        [S0, D, RESNORM] = diffusion_model_fit( ...
            Y, ...
            scheme, ...
            modelname = modelname,...
            fittingtechnique = fittingtechnique,...
            Nparam = 2,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );     


    case 'DKI'

        [S0, D, K, RESNORM] = diffusion_model_fit( ...
            Y, ...
            scheme, ...
            modelname = modelname,...
            fittingtechnique = fittingtechnique,...
            Nparam = 3,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );  


    case 'Sphere'

        [R, D, S0, RESNORM] = diffusion_model_fit( ...
            Y, ...
            scheme, ...
            modelname = modelname,...
            fittingtechnique = fittingtechnique,...
            Nparam = 3,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );


    case 'Ball+Sphere'

        [fs, R, Ds, Db, S0, RESNORM] = diffusion_model_fit( ...
            Y, ...
            scheme, ...
            modelname = modelname,...
            fittingtechnique = fittingtechnique,...
            Nparam = 5,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );

end


%% Save results

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting');

outf = fullfile(outputfolder, SampleName, modelname, schemename, fittingtechnique);
mkdir(outf)

% Meta data
if strcmp(fittingtechnique, 'LSQ')
    Meta = struct();
    Meta.beta0=beta0;
    Meta.lb=lb;
    Meta.ub=ub;
    Meta.lambda=lambda;
end



switch modelname

    case 'ADC'

        save(fullfile(outf, 'D.mat'), 'D');
        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM'); 
        save(fullfile(outf, 'Meta.mat'), 'Meta');  


    case 'DKI'

        save(fullfile(outf, 'D.mat'), 'D');
        save(fullfile(outf, 'K.mat'), 'K');
        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');    
        save(fullfile(outf, 'Meta.mat'), 'Meta');  

    case 'Sphere'

        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'D.mat'), 'D');    
        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');   
        save(fullfile(outf, 'Meta.mat'), 'Meta');    


    case 'Ball+Sphere'

        save(fullfile(outf, 'fs.mat'), 'fs');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'Ds.mat'), 'Ds');    
        save(fullfile(outf, 'Db.mat'), 'Db');  
        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');   
        save(fullfile(outf, 'Meta.mat'), 'Meta');    

end


% AIC
AIC = (nscheme)*log((RESNORM)/nscheme)-2*Nparam;
save(fullfile(outf, 'AIC.mat'), 'AIC'); 