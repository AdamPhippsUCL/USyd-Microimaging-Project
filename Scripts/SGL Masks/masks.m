% SCRIPT to auto generate tissue masks using 3DMGE, D, FA

clear;
projectfolder = pwd;

%% Sample and Image details

% Sample name
SampleName = '20250414_UQ6';
% SampleName = '20250407_UQ5';
SampleName = '20250224_UQ4';

% Use denoised data
UseDenoisedData = true;

% Imaging data folder 
switch UseDenoisedData
    case true
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN');     
    case false
        ImagingDataFolder = fullfile(projectfolder, 'Imaging Data', 'MAT');   
end

% 3D-MGE sequence
MGE_SeriesDescription = '3DMGE_20u';
MGE = load(fullfile(ImagingDataFolder, SampleName, MGE_SeriesDescription, 'avgImageArray.mat')).avgImageArray;

% dwFA
DTI_SeriesDescription = '40u_DtiSE_2012_SPOIL10% (20 micron)';
dwFA = load(fullfile(projectfolder, 'Outputs', 'Model Fitting', SampleName, 'DTI', DTI_SeriesDescription, 'dwFA.mat')).dwFA;


%%  Displaces images

switch SampleName
    case '20250407_UQ5'
        dx =-0;
        dy = -10;
        newMGE = zeros(size(MGE));
        newMGE(:,1:240+dx,1:640+dy)=MGE(:,1-dx:240,1-dy:640);
        MGE=newMGE;
    case '20250224_UQ4'
        disp('');
end

%% Generate masks

switch SampleName

    case '20250224_UQ4'

        MGElow = 5e-8;
        MGEhigh = 1.2e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);        
        

    case '20250407_UQ5'

        MGElow = 1e-7;
        MGEhigh = 3.4e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);


    case '20250414_UQ6'

        MGElow = 2e-8;
        MGEhigh = 1.1e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);
end


%% Display masks

displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = logical(GLANDULAR);
displaymasks(:,:,:,2) = logical(STROMA);
displaymasks(:,:,:,3) = logical(LUMEN);

sl=140;
cols = 400:630;
rows = 54:200;
figure
imshow(squeeze(MGE(sl,rows,cols)),[0 prctile(squeeze(MGE(sl,rows,cols)), 99.9, 'all')]);
figure
imshow(squeeze(dwFA(sl,rows,cols)),[0 prctile(squeeze(dwFA(sl,rows,cols)), 99.9, 'all')]);
figure
imshow(squeeze(MGE(sl,rows,cols)),[0 prctile(squeeze(MGE(sl,rows,cols)), 99.9, 'all')]);
hold on
mask = imshow(squeeze(displaymasks(sl,rows,cols,:)));
set(mask, 'AlphaData', 0.2)


%% Save masks

folder = fullfile(projectfolder, 'Outputs', 'Masks',SampleName, MGE_SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'LUMEN.mat'), 'LUMEN');
save(fullfile(folder, 'STROMA.mat'), 'STROMA');
save(fullfile(folder, 'GLANDULAR.mat'), 'GLANDULAR');