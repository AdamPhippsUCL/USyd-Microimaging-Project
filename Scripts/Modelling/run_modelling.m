% Script to run DW-MRI modelling on imaging data

clear;
projectfolder = pwd;

%% Image details

% Sample name
SampleName = '20250224_UQ4';

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
% modeltype = 'DKI';
modeltype = 'RDI - 1 compartment - 2 param';

% Scheme name
schemename = '20250224_UQ4 AllDELTA';
schemesfolder = fullfile(projectfolder, 'Schemes');
load(fullfile(schemesfolder, schemename));
nscheme = length(scheme);

% Correct for non-zero b0 value
effb0 = 10.3;
for indx = 1:nscheme
   if scheme(indx).bval==0
       continue
   end
   scheme(indx).bval = scheme(indx).bval-effb0;
end

% Fitting technique
fittingtechnique = 'LSQ';

% === MLP fitting

% Model folder
modelsfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models');

% === LSQ fitting

% Regularisation

R = 25; Rlb = 1; Rub = 100;

d = 1; dlb = 0.1; dub = 3;

dIC = 0.5; dIClb = 0.1; dICub = 3;

dEES = 1; dEESlb = 0.1; dEESub = 3;

fIC = 0.5; fIClb = 0; fICub = 1;

S0 = 1; S0lb = 0.8; S0ub = 1.2;

K = 0.5; Klb = 0; Kub = 3;

lambda0=1e-2;

switch modeltype

    case 'ADC'
        
        Nparam=2;
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

    case 'RDI - 2 compartment - 4 param'

        Nparam = 4;
        beta0 = [fIC,R,dIC,dEES];
        lb = [fIClb,Rlb,dIClb, dEESlb];
        ub = [fICub,Rub,dICub, dEESub];
        
    case 'RDKI - 2 compartment - 5 param'

        Nparam = 5;
        beta0 = [fIC,R,dIC,dEES, K];
        lb = [fIClb,Rlb,dIClb, dEESlb, Klb];
        ub = [fICub,Rub,dICub, dEESub, Kub];

    case 'RSI - 2 compartment - 4 param'

        Nparam = 4;
        beta0 = [0.4,d,0.6,d];
        lb = [fIClb,dlb,fIClb, dlb];
        ub = [fICub,dub,fICub, dub];
end

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
            Y(:,:,:,seriesindx) = (1/p(1))*((bimg./b0img)-p(2));
        
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

    case 'RDI - 1 compartment - 2 param'

        [R, dIC, RESNORM] = RDI_fit( ...
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


    case 'RDI - 2 compartment - 3 param'

        [fIC, R, d, RESNORM, fIC_stderr, R_stderr, d_stderr] = RDI_fit( ...
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
        

    case 'RDI - 2 compartment - 4 param'

        [fIC, R, dIC, dEES, RESNORM, fIC_stderr, R_stderr, dIC_stderr, dEES_stderr] = RDI_fit( ...
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


    case 'RDKI - 2 compartment - 5 param'

        [fIC, R, dIC, dEES, K, RESNORM, fIC_stderr, R_stderr, dIC_stderr, dEES_stderr, K_stderr] = RDI_fit( ...
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


    case 'RSI - 2 compartment - 4 param'

        [f1, d1, f2, d2, RESNORM, f1_stderr, d1_stderr, f2_stderr, d2_stderr] = RDI_fit( ...
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


    % case 'RDI - 3 compartment - fixed params'
    % 
    %     [fs, fg, fl] = RDI_fit( ...
    %         Y, ...
    %         scheme, ...
    %         modeltype = modeltype,...
    %         fittingtechnique = fittingtechnique,...
    %         modelfolder = modelfolder...
    %         );


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

        % switch averagedirections
        % 
        %     case true
        %         [ADC, S0] = calcADC(Y, [scheme(:).bval]);
        % 
        %     case false
        % 
        %         ADC = zeros([Ndirec, size(b0img)]);
        %         for direcindx = 1:Ndirec
        %             ADC(direcindx,:,:,:) = calcADC(squeeze(Y(direcindx,:,:,:,:)), [scheme(:).bval]);
        %         end
        % 
        % end


    % case 'No VASC VERDICT (AMICO)'
    % 
    %     [fIC, fEES, fVASC, R] = verdict_fit( ...
    %         Y, ...
    %         scheme, ...
    %         modeltype = modeltype,...
    %         fittingtechnique = fittingtechnique,...
    %         modelfolder = modelfolder...
    %         );
    % 
    % 
    % case {'MFP v1', 'MFP v2'}
    % 
    %     [R, D] = mfp_fit( ...
    %         Y, ...
    %         scheme,...
    %         modeltype = modeltype,...
    %         fittingtechnique=fittingtechnique,...
    %         modelfolder = modelfolder...
    %         );


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

        save(fullfile(outf, 'fIC_stderr.mat'), 'fIC_stderr');
        save(fullfile(outf, 'R_stderr.mat'), 'R_stderr');
        save(fullfile(outf, 'd_stderr.mat'), 'd_stderr');       
 
        % Cellularity
        C = fIC./(R.^3);
        C_stderr = C.*(fIC_stderr./fIC + 3*R_stderr./R);
        save(fullfile(outf, 'C.mat'), 'C');
        save(fullfile(outf, 'C_stderr.mat'), 'C_stderr');

    
    case 'RDI - 2 compartment - 4 param'

        save(fullfile(outf, 'Meta.mat'), 'Meta');   

        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');
        save(fullfile(outf, 'dEES.mat'), 'dEES');

        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');  

        save(fullfile(outf, 'fIC_stderr.mat'), 'fIC_stderr');
        save(fullfile(outf, 'R_stderr.mat'), 'R_stderr');
        save(fullfile(outf, 'dIC_stderr.mat'), 'dIC_stderr');
        save(fullfile(outf, 'dEES_stderr.mat'), 'dEES_stderr');

        % Cellularity
        C = fIC./(R.^3);
        C_stderr = C.*(fIC_stderr./fIC + 3*R_stderr./R);
        save(fullfile(outf, 'C.mat'), 'C');
        save(fullfile(outf, 'C_stderr.mat'), 'C_stderr');


    case 'RDKI - 2 compartment - 5 param'

        save(fullfile(outf, 'Meta.mat'), 'Meta');   

        save(fullfile(outf, 'fIC.mat'), 'fIC');
        save(fullfile(outf, 'R.mat'), 'R');
        save(fullfile(outf, 'dIC.mat'), 'dIC');
        save(fullfile(outf, 'dEES.mat'), 'dEES');
        save(fullfile(outf, 'K.mat'), 'K');

        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');  

        save(fullfile(outf, 'fIC_stderr.mat'), 'fIC_stderr');
        save(fullfile(outf, 'R_stderr.mat'), 'R_stderr');
        save(fullfile(outf, 'dIC_stderr.mat'), 'dIC_stderr');
        save(fullfile(outf, 'dEES_stderr.mat'), 'dEES_stderr');
        save(fullfile(outf, 'K_stderr.mat'), 'K_stderr');

        % Cellularity
        C = fIC./(R.^3);
        C_stderr = C.*(fIC_stderr./fIC + 3*R_stderr./R);
        save(fullfile(outf, 'C.mat'), 'C');
        save(fullfile(outf, 'C_stderr.mat'), 'C_stderr');

    case 'RDI - 3 compartment - fixed params'

        save(fullfile(outf, 'fs.mat'), 'fs');
        save(fullfile(outf, 'fg.mat'), 'fg');
        save(fullfile(outf, 'fl.mat'), 'fl');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');          
        save(fullfile(outf, 'Meta.mat'), 'Meta');    

    case 'ADC'

        save(fullfile(outf, 'ADC.mat'), 'ADC');
        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');    

    case 'DKI'

        save(fullfile(outf, 'S0.mat'), 'S0');
        save(fullfile(outf, 'D.mat'), 'D');
        save(fullfile(outf, 'K.mat'), 'K');
        save(fullfile(outf, 'RESNORM.mat'), 'RESNORM');            

end


% AIC
AIC = (nscheme)*log((RESNORM)/nscheme)-2*Nparam;
save(fullfile(outf, 'AIC.mat'), 'AIC'); 