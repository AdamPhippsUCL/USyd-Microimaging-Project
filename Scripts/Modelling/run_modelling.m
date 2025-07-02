% Script to run DW-MRI modelling on imaging data

clear;
projectfolder = pwd;

%% Image details

% Sample name
SampleName = '20250523_UQ8'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

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
modeltype = 'RDI - 2 compartment - 4 param (S0)'; % 'ADC', 'DKI', 'RDI - 2 compartment - 4 param (S0)'

% Scheme name
schemename = '20250224_UQ4 AllDELTA';
schemesfolder = fullfile(projectfolder, 'Schemes');
load(fullfile(schemesfolder, schemename));
nscheme = length(scheme);

% % Correct for non-zero b0 value
% effb0 = 10.3;
% for indx = 1:nscheme
%    if scheme(indx).bval==0
%        continue
%    end
%    scheme(indx).bval = scheme(indx).bval-effb0;
% end


% Fitting technique
fittingtechnique = 'LSQ';

% === MLP fitting

% Model folder
modelsfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models');

% === LSQ fitting

switch modeltype

    case 'ADC'
        
        Nparam=2;

        d=1; dlb = 0.2; dub = 2;
        S0=1; S0lb = 0.8; S0ub = 1.2;

        beta0=[S0,d];
        lb=[S0lb,dlb];
        ub=[S0ub,dub];

    case 'DKI'
        
        Nparam=3;
        beta0=[S0,d,K];
        lb=[S0lb,dlb,Klb];
        ub=[S0ub,dub,Kub];

    case 'RDI - 1 compartment - 2 param'

        Nparam = 2;
        beta0 = [R,d];
        lb = [Rlb,dlb];
        ub = [Rub,dub];

    case 'RDI - 2 compartment - 3 param'

        Nparam = 3;
        beta0 = [fIC,R,d];
        lb = [fIClb,Rlb,dlb];
        ub = [fICub,Rub,dub];

    case 'RDI - 2 compartment - 4 param (S0)'

        Nparam = 5;

        fIC = 0.5; fIClb = 0; fICub = 1;
        R = 6.5; Rlb = 6; Rub = 7;
        dIC = 0.5; dIClb = 0.1; dICub = 2;
        dEES = 1; dEESlb = 0.1; dEESub = 2;
        S0 = 1; S0lb = 0.8; S0ub = 1.2;

        beta0 = [fIC,R,dIC,dEES,S0];
        lb = [fIClb,Rlb,dIClb,dEESlb,S0lb];
        ub = [fICub,Rub,dICub,dEESub,S0ub];
        
end

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

% Define MLP model folder
modelfolder = fullfile(modelsfolder, modeltype, schemename);

switch modeltype

    case 'ADC'

        [S0, ADC, RESNORM] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );     


    case 'DKI'

        [S0, D, K, RESNORM] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );  


    case 'RDI - 1 compartment - 2 param (S0)'

        [R, dIC, S0, RESNORM] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );


    case 'RDI - 2 compartment - 3 param (S0)'

        [fIC, R, d, S0, RESNORM] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );
        

    case 'RDI - 2 compartment - 4 param (S0)'

        [fIC, R, dIC, dEES, S0, RESNORM] = RDI_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique = fittingtechnique,...
            modelfolder = modelfolder,...
            beta0=beta0,...
            lambda=lambda,...
            lb=lb,...
            ub=ub...
            );

   


end


%% Save results

outputfolder = fullfile(projectfolder, 'Outputs', 'Model Fitting');

outf = fullfile(outputfolder, SampleName, modeltype, schemename, fittingtechnique);
mkdir(outf)

% Meta data
if strcmp(fittingtechnique, 'LSQ')
    Meta = struct();
    Meta.beta0=beta0;
    Meta.lb=lb;
    Meta.ub=ub;
    Meta.lambda=lambda;
elseif strcmp(fittingtechnique, 'MLP')
    Meta = load(fullfile(modelfolder, 'Meta.mat')).Meta;
end


switch modeltype

    case 'ADC'

        save(fullfile(outf, 'ADC.mat'), 'ADC');
        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM'); 

    case 'RDI - 1 compartment - 2 param'

        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');    
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');   
        save(fullfile(outf, 'Meta.mat'), 'Meta');    


    case 'RDI - 2 compartment - 3 param'

        save(fullfile(outf, 'Meta.mat'), 'Meta'); 

        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'd.mat'), 'd');

        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');  

        % save(fullfile(outf, 'fIC_stderr.mat'), 'fIC_stderr');
        % save(fullfile(outf, 'R_stderr.mat'), 'R_stderr');
        % save(fullfile(outf, 'd_stderr.mat'), 'd_stderr');       
 
        % Cellularity
        C = fIC./(R.^3);
        save(fullfile(outf, 'C.mat'), 'C');
        % C_stderr = C.*(fIC_stderr./fIC + 3*R_stderr./R);
        % save(fullfile(outf, 'C_stderr.mat'), 'C_stderr');

    
    case 'RDI - 2 compartment - 4 param (S0)'

        save(fullfile(outf, 'Meta.mat'), 'Meta');   

        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');
        save(fullfile(outf, 'dEES.mat'), 'dEES');

        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');  

        % save(fullfile(outf, 'fIC_stderr.mat'), 'fIC_stderr');
        % save(fullfile(outf, 'R_stderr.mat'), 'R_stderr');
        % save(fullfile(outf, 'dIC_stderr.mat'), 'dIC_stderr');
        % save(fullfile(outf, 'dEES_stderr.mat'), 'dEES_stderr');

        % Cellularity
        C = fIC./(R.^3);
        save(fullfile(outf, 'C.mat'), 'C');
        % C_stderr = C.*(fIC_stderr./fIC + 3*R_stderr./R);
        % save(fullfile(outf, 'C_stderr.mat'), 'C_stderr');


end


% AIC
AIC = (nscheme)*log((RESNORM)/nscheme)-2*Nparam;
save(fullfile(outf, 'AIC.mat'), 'AIC'); 