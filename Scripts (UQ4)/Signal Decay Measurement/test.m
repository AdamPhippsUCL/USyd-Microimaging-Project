% TEST script to measure R and D for each compartment (E, S, L)

clear;

% Project folder
projectfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project";

samplename = '20250224_UQ4';

%% Base image (3DMGE)

seriesdescription = '3DMGE_20u';

% Load image
baseImg = load(fullfile(projectfolder, 'Imaging Data', 'MAT DN', samplename, seriesdescription, 'avgImageArray.mat')).avgImageArray;
szbase = size(baseImg);


%% MASKS

maskfolder = fullfile(projectfolder, 'Outputs', 'Masks', samplename, seriesdescription);
GLANDULAR = load(fullfile(maskfolder, 'GLANDULAR.mat')).GLANDULAR;
STROMA = load(fullfile(maskfolder, 'STROMA.mat')).STROMA;
LUMEN = load(fullfile(maskfolder, 'LUMEN.mat')).LUMEN;


%% Display masks

displaymasks = zeros([size(GLANDULAR), 3]);
displaymasks(:,:,:,1) = logical(GLANDULAR);
displaymasks(:,:,:,2) = logical(STROMA);
displaymasks(:,:,:,3) = logical(LUMEN);

sl = 128;
figure
imshow(squeeze(baseImg(sl,:,:)),[]);
hold on
mask = imshow(squeeze(displaymasks(sl,:,:,:)));
set(mask, 'AlphaData', 0.2)



%% Make sample mask

switch samplename

    case '20250224_UQ4'

        % Cylinder centred at (128, 114)
        
        samplemask = zeros(size(baseImg));
        
        [Xs, Ys] = meshgrid( ...
           1:size(baseImg,2), ...
           1:size(baseImg,1) ...
            );
        
        samplemask(:,:,:) = repmat((Xs-128).^2 + (Ys-114).^2 <75^2, 1, 1, size(baseImg,3));

end

%% LOAD LOW RES DW-MRI images

schemename = '20250224_UQ4 AllDELTA';
schemesfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\DW-MRI-Modelling\Schemes";
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


UseDenoisedData = true;

SeriesDescriptions = {
    'SE_b0_SPOIL5% (640 micron)',...
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
};

Nimg = length(SeriesDescriptions)-1;


%% Data preprocessing

% == Load images and dinfo

ImageArrays = struct();
DINFOS = struct();

for seriesindx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{seriesindx};

    % Load image and dinfo
    switch UseDenoisedData
        case true
            thisfolder = fullfile(projectfolder, 'Imaging Data', 'MAT DN', samplename, SeriesDescription);
        case false
            thisfolder = fullfile(projectfolder, 'Imaging Data', 'MAT', samplename, SeriesDescription);
    end
    ImageArray = load(fullfile(thisfolder, 'axialImageArray.mat')).ImageArray;
    dinfo = load(fullfile(thisfolder, 'axialdinfo.mat')).dinfo;

    % Append to structures
    DINFOS(seriesindx).dinfo = dinfo;
    ImageArrays(seriesindx).ImageArray = ImageArray;

end


szimg = size(ImageArray, 1:3);
IMGS = zeros([Nimg, szimg]);

% seriesindx = 1 used for data normalisation!
img = ImageArrays(1).ImageArray;
dinfo = DINFOS(1).dinfo;
bvals = [dinfo(:).DiffusionBValue];
b0bools = (bvals==0);
b0imgs = img(:,:,:,b0bools);
b0img = mean(b0imgs,4);
meanb0 = mean(b0img(:));


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
    IMGS(seriesindx-1,:,:,:) = (meanb0/thismeanb0)*(bimg./b0img); 

end




%% Construct composition image

% Choose low resolution image size
szmap = szimg;
COMPOSITION = zeros([szmap, 3]);
ResFactor = szbase./szmap;
Nvoxel = prod(ResFactor);

for rindx = 1:szmap(1)
    for cindx = 1:szmap(2)
        for slindx = 1:szmap(3)

            baserows = ((rindx-1)*ResFactor(1)+1:rindx*ResFactor(1));
            basecols = ((cindx-1)*ResFactor(2)+1:cindx*ResFactor(2));
            baseslices = ((slindx-1)*ResFactor(3)+1:slindx*ResFactor(3));

            % Test in sample
            if ~all(logical(samplemask(baserows, basecols, baseslices)), "all")
                continue
            end

            % COMPOSITION         
            COMPOSITION(rindx, cindx, slindx, 1) = sum(double(STROMA(baserows, basecols, baseslices)), "all")/Nvoxel;
            COMPOSITION(rindx, cindx, slindx, 2) = sum(double(GLANDULAR(baserows, basecols, baseslices)), "all")/Nvoxel;
            COMPOSITION(rindx, cindx, slindx, 3) = sum(double(LUMEN(baserows, basecols, baseslices)), "all")/Nvoxel;        

        end
    end
end

STROMAvals = squeeze(COMPOSITION(:,:,:,1));
GLANDULARvals = squeeze(COMPOSITION(:,:,:,2));
LUMENvals = squeeze(COMPOSITION(:,:,:,3));


%% Measurements

% Flatten COMPOSITION and remove zeros
composition = reshape(COMPOSITION,[prod(szmap), 3]);
bool =  sum(composition,2)>0;
% bool = (composition(:,3)>0.1);
composition = composition(bool, :);

% Flatten IMGS and remove zeros (from COMPOSITION)
imgs = reshape(IMGS, [Nimg, prod(szimg)]);
imgs = imgs(:,bool);



%% Test example 

% indxs = 1:340;
% 
% y = imgs(:,indxs);
% y = reshape(y, [numel(y),1]);
% 
% X = zeros(size(y,1), 6);
% 
% for ischeme = 1:nscheme/2
%     X(ischeme:nscheme/2:end,1)= scheme(2*ischeme).bval;
%     X(ischeme:nscheme/2:end,2)= scheme(2*ischeme).delta;
%     X(ischeme:nscheme/2:end,3)= scheme(2*ischeme).DELTA;
% end
% 
% comp = composition(indxs, :);
% for i = 1:length(indxs)
%     X(nscheme/2*(i-1)+1:nscheme/2*i,4:6)=repmat(composition(indxs(i),:),nscheme/2,1);
% end
% 
% func = @(b,X) test_func(b,X, comp);
% beta0 = [10,2,10,2,35,2];
% lower_bounds = [0.1,0.1,0.1,0.1,0.1,0.1];
% upper_bounds = [40,3,40,3,40,3];
% [beta, resnorm] = lsqcurvefit(@test_func, beta0, X, y, lower_bounds, upper_bounds);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% beta = lsqcurvefit(func, beta0, X, y, lower_bounds, upper_bounds);

%% New test...

indxs = 1:size(composition, 1);

% % Bootstrap samples
% n = length(indxs);
% Nboot = 1001;
% bootindxs = randi(n, n, Nboot);

% Optimisation
func = @(b, X) test_func2(b,X);
beta0 = [0.7,0.5,0.1];
lb = [0.0,0.0,0.0];
ub = [1,1,1];

% Lumen diffusivity
Dl = 0.0020;
err = 0.001;
    
% signals = zeros(3, Nimg,Nboot); % BOOTSTRAPPING
signals = zeros(3, Nimg, 2); % STANDARD ERROR
for imgindx = 1:Nimg
      
    X = composition(indxs,:);
    y = transpose(imgs(imgindx, indxs));

    % Set lumen prior
    Sl = exp(-scheme(2*imgindx).bval*Dl);
    lb(3)=Sl-err;
    ub(3)=Sl+err;


    % BOOTSTRAPPING METHOD

    % bootsignals = zeros(Nboot, 3);
    % 
    % for bootindx = 1:Nboot
    % 
    %     thisX = X(bootindxs(:,bootindx),:);
    %     thisy = y(bootindxs(:,bootindx),:);
    % 
    %     [beta_fit, resnorm, residual, exitflag, output, lambda, jacobian]=lsqcurvefit(func, beta0, thisX, thisy, lb, ub);
    %     bootsignals(bootindx,:) = beta_fit;
    % 
    % end
    % 
    % 
    % % signals(:,imgindx,1) = mean(bootsignals, 1);
    % % signals(:,imgindx,2) = prctile(bootsignals,5,1);
    % % signals(:,imgindx,3) = prctile(bootsignals,95,1);
    % 
    % signals(1,imgindx,:) = bootsignals(:,1);
    % signals(2,imgindx,:) = bootsignals(:,2);
    % signals(3,imgindx,:) = bootsignals(:,3);

    % STANDARD ERROR
    [beta_fit, resnorm, residual, exitflag, output, lambda, jacobian]=lsqcurvefit(func, beta0, X, y, lb, ub);

    % MSE
    mse = resnorm / (length(y) - length(beta0));

    % Covariance matrix
    cov_matrix = mse * inv(jacobian' * jacobian);
    
    % Standard errors (square root of diagonal elements)
    stderr = sqrt(diag(cov_matrix));

    signals(:,imgindx,1) = beta_fit;
    signals(:,imgindx,2) = stderr;

    % 
    % y_pred = func(beta_fit, X);
    % figure
    % scatter(y_pred, y)


end

%% Display

% % DISPLAY RESULTS FROM BOOTSTRAPPING
% 
% bvals = [scheme(2:2:end).bval];
% bshift = 25;
% width = 15;
% figure
% h1=boxplot(transpose(squeeze(signals(1,:,:))), 'Positions', bvals-1*bshift, 'Colors', 	[0.4660 0.6740 0.1880], 'Widths', width, 'Symbol', 'k*');
% hold on
% h2=boxplot(transpose(squeeze(signals(2,:,:))), 'Positions', bvals-0*bshift, 'Colors', [0.8500 0.3250 0.0980], 'Widths', width, 'Symbol', 'k*');
% h3=boxplot(transpose(squeeze(signals(3,:,:))), 'Positions', bvals+1*bshift, 'Colors', [0 0.4470 0.7410], 'Widths', width, 'Symbol', 'k*');
% xticks(bvals); 
% xticklabels(bvals)
% ylim([-4,0])
% xlabel('b-value')
% ylabel('Predicted normalized signal')
% 
% % Find the outlier markers
% outliers = findobj(gca, 'Tag', 'Outliers');
% % Reduce the size and linewidth of the outliers
% for i = 1:length(outliers)
%     outliers(i).MarkerSize = 1;   % Reduce size (default ~6)
%     outliers(i).LineWidth = 0.2;  % Reduce linewidth (default ~0.75-1)
% end
% 
% 
% plot(bvals-1*bshift, median(squeeze(signals(1,:,:)), 2), '--', color=[0.4660 0.6740 0.1880], HandleVisibility='off');
% plot(bvals-0*bshift, median(squeeze(signals(2,:,:)), 2), '--', color=[0.8500 0.3250 0.0980], HandleVisibility='off');
% plot(bvals+1*bshift, median(squeeze(signals(3,:,:)), 2), '--', color=[0 0.4470 0.7410], HandleVisibility='off');
% 
% 
% % Get handles for boxes (for legend)
% boxHandles1 = findobj(h1, 'Tag', 'Box'); medianHandles1 = findobj(gca, 'Tag', 'Median');
% set(boxHandles1, 'LineWidth', 1); set(medianHandles1, 'LineWidth', 1);
% boxHandles2 = findobj(h2, 'Tag', 'Box'); medianHandles2 = findobj(gca, 'Tag', 'Median');
% set(boxHandles2, 'LineWidth', 1); set(medianHandles2, 'LineWidth', 1);
% boxHandles3 = findobj(h3, 'Tag', 'Box'); medianHandles3 = findobj(gca, 'Tag', 'Median');
% set(boxHandles3, 'LineWidth', 1); set(medianHandles3, 'LineWidth', 1);
% 
% % Create legend handles (use a dummy line for the legend)
% l1 = plot(NaN, NaN, 's', 'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor', [0.4660 0.6740 0.1880]); % Dummy handle for blue group
% l2 = plot(NaN, NaN, 's', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', [0.8500 0.3250 0.0980]); % Dummy handle for red group
% l3 = plot(NaN, NaN, 's', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410]); % Dummy handle for green group
% 
% % Add a legend
% legend([l1, l2, l3], {'Stroma', 'Glandular', 'Lumen'}, 'Location', 'northeast');


% DISPLAY STANDARD ERROR RESULTS

figure



% Deltas = [scheme(indices).DELTA];
% % [Ds, I] = sort(Deltas);
% Ds=Deltas;
% I=1:length(Ds);
% Dshift=2;

% Short Delta
indices = 2:2:10;
s = signals(:,1:5,:);
bvals = [scheme(indices).bval];
bshift = 10;

errorbar(bvals-1*bshift, s(1,:,1), s(1,:,2), '--*', color=[0.4660 0.6740 0.1880], DisplayName='Stroma (Short Delta)');
hold on
errorbar(bvals-1*bshift, s(2,:,1), s(2,:,2), '--*', color=[0.8500 0.3250 0.0980], DisplayName='Glandular (Short Delta)');
errorbar(bvals-1*bshift, s(3,:,1), s(3,:,2), '--*', color=[0 0.4470 0.7410], DisplayName = 'Lumen (Short Delta)');


% Long Delta
indices = 12:2:20;
s = signals(:,6:10,:);

errorbar(bvals+1*bshift, s(1,:,1), s(1,:,2), '--.', color=[0.4660 0.6740 0.1880], DisplayName='Stroma (Long Delta)');
hold on
errorbar(bvals+1*bshift, s(2,:,1), s(2,:,2), '--.', color=[0.8500 0.3250 0.0980], DisplayName='Glandular (Long Delta)');
errorbar(bvals+1*bshift, s(3,:,1), s(3,:,2), '--.', color=[0 0.4470 0.7410], DisplayName = 'Lumen (Long Delta)');



xticks(bvals); 
xticklabels(bvals)
ylim([-0.1,0.8])
xlim([800,2200])
xlabel('b-value')
grid on
% errorbar(Ds-1*Dshift, s(1,I,1), s(1,I,2), '--*', color=[0.4660 0.6740 0.1880], DisplayName='Stroma');
% hold on
% errorbar(Ds+0*Dshift, s(2,I,1), s(2,I,2), '--*', color=[0.8500 0.3250 0.0980], DisplayName='Glandular');
% errorbar(Ds+1*Dshift, s(3,I,1), s(3,I,2), '--*', color=[0 0.4470 0.7410], DisplayName = 'Lumen');
% legend;
% xlabel('Delta')
ylabel('Predicted normalized signal')
% title(T)

legend;

% % Create legend handles (use a dummy line for the legend)
% l1 = plot(NaN, NaN, 's', 'MarkerFaceColor', [0.4660 0.6740 0.1880], 'MarkerEdgeColor', [0.4660 0.6740 0.1880]); % Dummy handle for blue group
% l2 = plot(NaN, NaN, 's', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerEdgeColor', [0.8500 0.3250 0.0980]); % Dummy handle for red group
% l3 = plot(NaN, NaN, 's', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', [0 0.4470 0.7410]); % Dummy handle for green group
% 
% % Add a legend
% legend([l1, l2, l3], {'Stroma', 'Glandular', 'Lumen'}, 'Location', 'northeast');

ax=gca();
ax.FontSize=14;


%% INITIALISE RESULTS STRUCTURE

RESULTS = struct();


%% MODEL FITTING

% delete 'results_tab.xlsx'

component = 'G';


switch component
    case 'S'
        S = squeeze(signals(1,:,1));
    case 'G'
        S = squeeze(signals(2,:,1));
    case 'L'
        S = squeeze(signals(3,:,1));
end

function Y_out = RDI_model_func(b,X, modeltype, RDI_model)   

    N = size(X,1);
    Y_out = zeros([N,1]);

    for i=1:N
        Y_out(i) = RDI_model(b,X(i,:),modeltype=modeltype);
    end

end

% Model fitting
modeltype = 'RDI - 2 compartment - 4 param';
% modeltype = 'ADC';

R = 30;
Rlb = 1;
Rub = 100;

d = 1;
dlb = 0.1;
dub = 3;

fIC = 0.5;
fIClb = 0;
fICub = 1;

S0 = 1;
S0lb = 0.8;
S0ub = 1.2;

K = 1;
Klb = 0;
Kub = 3;

lambda0=5e-2;

switch modeltype

    case 'ADC'

        Nparam = 2;
        beta0 = [S0, d];
        lb = [S0lb, dlb];
        ub = [S0ub, dub];

    case 'DKI'

        Nparam = 3;
        beta0 = [S0, d, K];
        lb = [S0lb, dlb, Klb];
        ub = [S0ub, dub, Kub];

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
        beta0 = [fIC,R,d,d];
        lb = [fIClb,Rlb,dlb, dlb];
        ub = [fICub,Rub,dub, dub];

    case 'RDKI - 2 compartment - 5 param'

        Nparam = 5;
        beta0 = [fIC,R,d,d, K];
        lb = [fIClb,Rlb,dlb, dlb, Klb];
        ub = [fICub,Rub,dub, dub, Kub];


    case 'RSI - 2 compartment - 4 param'

        Nparam = 4;
        beta0 = [0.4,d,0.6,d];
        lb = [fIClb,dlb,fIClb, dlb];
        ub = [fICub,dub,fICub, dub];

end


lambda = lambda0*ones(1,Nparam);

Y = ones([1,1,1,nscheme]);
Y(:,:,:,2:2:end) = S;

schemename = schemename;
fittingtechnique = 'LSQ';



modelfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models', modeltype, schemename);

clear outputs;
clear stderr;
[outputs{1:Nparam}, resnorm, stderr{1:Nparam}] = RDI_fit( ...
    Y, ...
    scheme, ...
    modeltype=modeltype, ...
    modelfolder=modelfolder, ...
    fittingtechnique=fittingtechnique,...
    beta0=beta0,...
    lb=lb,...
    ub=ub,...
    lambda=lambda);

model_params = cell2mat(outputs);

% Predictions
Y_pred = zeros(size(Y));
for indx = 1:length(Y)
    bval = scheme(indx).bval;
    delta = scheme(indx).delta;
    DELTA = scheme(indx).DELTA;
    Y_pred(indx) = RDI_model( model_params, [bval,delta,DELTA], modeltype = modeltype);
end

% Residual error
reserror = sum((Y_pred(2:2:end)-Y(2:2:end)).^2);

% AIC
N=nscheme/2;
AIC = (N)*log(reserror/N)+2*Nparam;

figure
scatter(squeeze(Y), squeeze(Y_pred))
hold on
plot(squeeze(Y), squeeze(Y))
title(modeltype)

n = length(RESULTS)+1;
if ~numel(fieldnames(RESULTS))
    n = 1;
end
RESULTS(n).Component = component;
RESULTS(n).ModelType = modeltype;
RESULTS(n).ModelParams = model_params;
RESULTS(n).StdErr = cell2mat(stderr);
RESULTS(n).ResidualError = reserror;
RESULTS(n).AIC=AIC;

tab=struct2table(RESULTS);
disp(tab);
% writetable(tab, 'results_tab.xlsx')


% % == 2 compartment
% 
% % Model fitting
% modeltype = 'RDI - 2 compartment - 4 param';
% Nparam = 4;
% schemename = schemename;
% 
% % MLP fitting
% 
% Y = ones([1,1,1,nscheme]);
% Y(:,:,:,2:2:end) = S;
% 
% fittingtechnique = 'MLP';
% modelfolder = fullfile(projectfolder, 'Scripts', 'RDI', fittingtechnique, 'models', modeltype, schemename);
% [fIC, fEES, R, dIC, dEES]= RDI_fit(Y, scheme, modeltype=modeltype, modelfolder=modelfolder, fittingtechnique=fittingtechnique)
% 
% 
% % LSQ FITTING
% 
% Y = ones([nscheme,1]);
% Y(2:2:end) = S;
% 
% RDI2_func = @(b,X) RDI_model_func(b, X, modeltype, @RDI_model);
% 
% X = zeros(nscheme, 3);
% X(:,1) = [scheme(:).bval];
% X(:,2) = [scheme(:).delta];
% X(:,3) = [scheme(:).DELTA];
% 
% beta0 = [0.5, 20, 2, 2];
% lb = [0, 10, 0.1, 0.1];
% ub = [1, 20,1, 1];
% 
% beta = lsqcurvefit(RDI2_func, beta0, X, Y, lb, ub);
% fIC = beta(1)
% R = beta(2)
% dIC = beta(3)
% dEES = beta(4)
% 
% 
% % Predictions
% Y_pred = zeros(size(Y));
% for indx = 1:length(Y)
%     bval = scheme(indx).bval;
%     delta = scheme(indx).delta;
%     DELTA = scheme(indx).DELTA;
%     Y_pred(indx) = RDI_model( [fIC, R, dIC, dEES], [bval,delta,DELTA], modeltype = modeltype);
% end
% 
% % Residual error
% reserror = sum((Y_pred(2:2:end)-Y(2:2:end)).^2);
% 
% % AIC
% RDI2_AIC = (nscheme/2)*log(2*reserror/nscheme)+2*Nparam
% 
% 
% figure
% scatter(squeeze(Y), squeeze(Y_pred))
% hold on
% plot(squeeze(Y), squeeze(Y))
% title('RDI - 2 compartment - 4 param')


%% DKI MODEL FITTING

% S = squeeze(signals(1,:,1));
% 
% Y = ones([nscheme,1]);
% Y(2:2:end) = S;
% 
% Nparam = 3;
% 
% 
% function Y_out = DKI_model_func(b,X, DKI_model)   
%     N = size(X,1);
%     Y_out = zeros([N,1]);
% 
%     for i=1:N
%         Y_out(i) = DKI_model(b,X(i,:));
%     end
% 
% end
% 
% 
% DKI_func = @(b,X) DKI_model_func(b, X, @DKI_model);
% 
% X = zeros(nscheme, 1);
% X(:,1) = [scheme(:).bval];
% 
% beta0 = [1, 0.002, 0.5];
% lb = [0.8,0.0001, 0];
% ub = [1.2, 0.1, 3];
% 
% beta = lsqcurvefit(DKI_func, beta0, X, Y, lb, ub);
% S0 = beta(1)
% D = beta(2)
% K = beta(3)
% 
% 
% % Predictions
% Y_pred = zeros(size(Y));
% for indx = 1:length(Y)
%     bval = scheme(indx).bval;
%     Y_pred(indx) = DKI_model( [S0, D, K], [bval]);
% end
% 
% figure
% scatter(squeeze(Y), squeeze(Y_pred))
% hold on
% plot(squeeze(Y), squeeze(Y))
% title('DKI')
% 
% 
% DKI_AIC = (nscheme/2)*log(2*reserror/nscheme)+2*Nparam