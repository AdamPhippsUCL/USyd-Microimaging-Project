% Script to calculate AIC of each model on signals from individual samples

clear;
projectfolder = pwd;

%% Load modelling results

RESULTS = load(fullfile(projectfolder, 'Outputs', 'Signal Measurement', 'RESULTS.mat')).RESULTS;


%% Sample and scheme details

outputfolder = fullfile(projectfolder, 'Outputs', 'Signal Measurement'); 

samplename = '20250224_UQ4';
schemename = '20250224_UQ4 AllDELTA';

% Measured signals and scheme
signals = load(fullfile(outputfolder, samplename, schemename, 'signals.mat')).signals;
scheme = load(fullfile(outputfolder, samplename, schemename, 'scheme.mat')).scheme;
nscheme = length(scheme);


%% Calculate AIC for each model entry

for rindx = 1:length(RESULTS)

   % Tissue component
   tissuecomponent = RESULTS(rindx).Component;

   switch tissuecomponent
       case 'S'
           compindx = 1;
       case 'G'
           compindx = 2;
       case 'L'
           compindx = 3;
   end

   % model type
   modeltype = RESULTS(rindx).ModelType;

   % parameters 
   params = RESULTS(rindx).ModelParams;
   Nparam = length(params);


   % Simulate signals using parameters
   simsignals = zeros(1, nscheme);
   for ischeme = 1:nscheme
        
       bval = scheme(ischeme).bval;
       delta = scheme(ischeme).delta;
       DELTA = scheme(ischeme).DELTA;

       simsignals(ischeme) = RDI_model(params, [bval, delta, DELTA], modeltype=modeltype);


   end
   
   f=figure;
   plot([scheme(:).bval], signals(compindx,:,1), '*');
   hold on
   plot([scheme(:).bval], simsignals, '*');
   close(f)

   % Residual error (sum of squares)
   resnorm = sum( (signals(compindx,:,1)-simsignals).^2);


   % AIC
   AIC = nscheme*log(resnorm/nscheme) + 2*Nparam;

   % Append to RESULTS
   RESULTS(rindx).(['AIC_' samplename(end-2:end)]) = AIC;
    

end