% Script to define and make scheme

% Scheme: [delta,
%          DELTA, 
%          b, 
%          TE, 
%          TR, 
%          NSA,  (Number of signal averages for b=0 image)
%          Rav   (Additional signal averaging factor for b>0 image)
%           ]  

%% Define Scheme

% UQ 1
schemename = 'UQ Scheme v2';

V01 = [1, 2, 0, 26, 250, 4, 1.5];
V1 = [1, 20, 1000, 26, 250, 4, 1.5];
V02 = [1, 2, 0, 23.5, 250, 4, 1.5];
V2 = [1, 17.5, 1500, 23.5, 250, 4, 1.5];
V03 = [1, 2, 0, 21, 250, 4, 1.5];
V3 = [1, 15, 2000, 21, 250, 4, 1.5];
V04 = [1, 2, 0, 18.5, 250, 4, 1.5];
V4 = [1, 12.5, 2500, 18.5, 250, 4, 1.5];
V05 = [1, 2, 0, 16, 250, 4, 1.5];
V5 = [1, 10, 3000, 16, 250, 4, 1.5];

Vs = [...
    V01; V1;...
    V02; V2;...
    V03; V3;...
    V04; V4;...
    V05; V5;...            
    ];




%% Build and save scheme

scheme = BuildScheme(Vs, schemename);

% Save scheme
schemefolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Scripts\VERDICT\Schemes";
save([char(schemefolder) '\' schemename '.mat'], 'scheme')
