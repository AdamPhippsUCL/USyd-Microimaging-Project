% Script to define and make DW-MRI scheme structure

clear;
projectfolder = pwd;

%% Define Scheme

schemename = '20250224_UQ4 AllDELTA';

% Scheme: [delta, DELTA, b]  

V0 = [1, 2, 0];
V1 = [2, 15, 1006];
V2 = [2, 20, 1259];
V3 = [2, 30, 1514];
V4 = [2, 40, 1769];
V5 = [2, 50, 2023];
V6 = [2, 40, 1015];
V7 = [2, 60, 1271];
V8 = [2, 80, 1528];
V9 = [2, 100, 1785];
V10 = [2, 120, 2042];

Vs = [V0; V1; V2; V3; V4; V5; V6; V7; V8; V9; V10];


%% Build and save scheme

scheme = BuildScheme(Vs, schemename);

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
