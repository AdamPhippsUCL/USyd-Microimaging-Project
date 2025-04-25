% Script to define and make DW-MRI scheme structure

clear;
projectfolder = pwd;

%% Define Scheme

schemename = '20250224_UQ4 AllDELTA';

% Scheme: [delta, DELTA, b]  (b = mean(eff b values) - effSTEb0) 



V0 = [1, 2, 0];
V1 = [2, 15, 1000];
V2 = [2, 20, 1250];
V3 = [2, 30, 1500];
V4 = [2, 40, 1750];
V5 = [2, 50, 2000];
V6 = [2, 40, 1000];
V7 = [2, 60, 1250];
V8 = [2, 80, 1500];
V9 = [2, 100, 1750];
V10 = [2, 120, 2000];

% V0 = [1, 2, 0];
% V1 = [2, 15, 997];
% V2 = [2, 20, 1246];
% V3 = [2, 30, 1492];
% V4 = [2, 40, 1737];
% V5 = [2, 50, 1983];
% V6 = [2, 40, 983];
% V7 = [2, 60, 1222];
% V8 = [2, 80, 1461];
% V9 = [2, 100, 1700];
% V10 = [2, 120, 1938];

Vs = [V0; V1; V2; V3; V4; V5; V6; V7; V8; V9; V10];

effb0vals = [10.3, 9.05, 13.5, 22.5, 31.5, 40.4, 31.4, 49.4,67.3, 85.2,103.1];


%% Build and save scheme

scheme = BuildScheme(Vs, schemename);


effb0vals = num2cell(effb0vals);
[scheme(:).effb0vals] = deal(effb0vals{:});

% Save scheme
schemesfolder = fullfile(projectfolder, 'Schemes');
save(fullfile(schemesfolder, [schemename '.mat']));

% Build a scheme structure
function scheme = BuildScheme(scans, schemename)

for scanIndx = 1:size(scans,1)

    delta = scans(scanIndx, 1);
    Delta = scans(scanIndx,2);
    bval = scans(scanIndx,3);

    scheme(scanIndx).delta = delta;
    scheme(scanIndx).DELTA = Delta ;
    scheme(scanIndx).bval = bval;
    scheme(scanIndx).G = stejskal(delta,Delta,bval=bval);
    scheme(scanIndx).schemename = schemename;

end
    
end
