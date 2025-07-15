% SCRIPT to auto generate tissue masks using 3DMGE, D, FA

clear;
projectfolder = pwd;

%% Sample and Image details

% Sample name
SampleName = '20250414_UQ6'; % '20250224_UQ4', '20250407_UQ5', '20250414_UQ6', '20250522_UQ7', '20250523_UQ8', '20250524_UQ9'

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
        MGEhigh = 1.20e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);        
        

    case '20250407_UQ5'

        MGElow = 0.5e-7;
        MGEhigh = 3.5e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);


    case '20250414_UQ6'

        MGElow = 2.5e-8;
        MGEhigh = 1.05e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);

    case '20250522_UQ7'

        MGElow = 1.0e-8;
        MGEhigh = 4.5e-8;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);


    case '20250523_UQ8'

        MGElow = 4e-8;
        MGEhigh = 1.6e-7;
        dwFAlow = 14e-5;

        STROMA = (dwFA>dwFAlow).*and(MGE<MGEhigh, MGE>MGElow);
        LUMEN = (MGE>MGEhigh);
        GLANDULAR = and(~logical(STROMA), ~logical(LUMEN)).*(MGE>MGElow);


    case '20250524_UQ9'

        MGElow = 3.0e-8;
        MGEhigh = 1.05e-7;
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

sl=120;
cols = 1:640;%20:620;
rows = 1:240;%35:210;
figure
imshow(squeeze(MGE(sl,rows,cols)),[0 prctile(squeeze(MGE(sl,rows,cols)), 99.9, 'all')]);
figure
imshow(squeeze(dwFA(sl,rows,cols)),[0 prctile(squeeze(dwFA(sl,rows,cols)), 99.9, 'all')]);
figure
imshow(squeeze(MGE(sl,rows,cols)),[0 prctile(squeeze(MGE(sl,rows,cols)), 99.9, 'all')]);
hold on
mask = imshow(squeeze(displaymasks(sl,rows,cols,:)));
set(mask, 'AlphaData', 0.2)

% Axial
sl=74;
figure
imshow(squeeze(MGE(:,:,sl)),[0 prctile(squeeze(MGE(:,:,sl)), 99.9, 'all')]);
hold on
mask = imshow(squeeze(displaymasks(:,:,sl,:)));
set(mask, 'AlphaData', 0.2)

%% Save masks

folder = fullfile(projectfolder, 'Outputs', 'Masks',SampleName, MGE_SeriesDescription);
mkdir(folder);
save(fullfile(folder, 'LUMEN.mat'), 'LUMEN');
save(fullfile(folder, 'STROMA.mat'), 'STROMA');
save(fullfile(folder, 'GLANDULAR.mat'), 'GLANDULAR');