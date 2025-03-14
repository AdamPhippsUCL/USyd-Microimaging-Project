% Script to save DTI sequence information as structure

%% Initial details

% Imaging data folder
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Image type
ImageType = 'MAT';

% Sample name
SampleName = '20250224_UQ4';

% Series Description
SeriesDescriptions = {
    'SE_b0_SPOIL10%',...
    'STEAM_DELTA_15',...
    'STEAM_DELTA_20',...
    'STEAM_DELTA_30',...
    'STEAM_DELTA_40',...
    'STEAM_DELTA_50',...
    '40u_DtiSTE_SPOIL5%'
};


%% Define structure

for indx = 1:length(SeriesDescriptions)

    SeriesDescription = SeriesDescriptions{indx};

    DTIstruct = struct();
    
    switch SampleName
    
        case '20241128_UQ1'
    
            switch SeriesDescription

                case '40u_verdict_seq1_v2'
    
                    DiffusionBValue = [0, 0, 0, 0,...
                                    1000, 1000, 1000,...
                                    1000, 1000, 1000];
    
                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 40.874806178439, 40.874806178439,...
                        1018.51474328831, 1003.52893354527, 1000.48691656872, 1025.2057821385,... 
                        1006.216319282, 1012.68675629359...                       
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  


                case '40u_verdict_seq2_v2'
    
                    DiffusionBValue = [0, 0, 0, 0,...
                                    1500, 1500, 1500,...
                                    1500, 1500, 1500];
    
                    DiffusionEffBValue = [...
                            40.874806178439, 40.874806178439, 40.874806178439, 40.874806178439,... 
                            1512.53599635572, 1500.86509101271, 1500.63692461532, 1522.21932326803,...
                            1500.60635346786, 1509.76052269534...                  
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  


                case '40u_verdict_seq3_v2'
    
                    DiffusionBValue = [0, 0, 0, 0,...
                                    2000, 2000, 2000,...
                                    2000, 2000, 2000];
    
                    DiffusionEffBValue = [...
                            40.874806178439, 40.874806178439, 40.874806178439, 40.874806178439,... 
                            2010.91724442071, 2001.07945135017, 2000.79448041395, 2017.90004704977,... 
                            2000.75629824247, 2005.59766781398...                
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];  


                case '40u_verdict_seq4_v2'
    
                    DiffusionBValue = [0, 0, 0, 0,...
                                    2500, 2500, 2500,...
                                    2500, 2500, 2500];
    
                    DiffusionEffBValue = [...
                            40.874806178439, 40.874806178439, 40.874806178439, 40.874806178439,... 
                            2508.40464778432, 2501.3237220966, 2500.97402049583, 2511.50113943203,...
                            2500.92716530889, 2499.64680644115...              
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 


                case '40u_verdict_seq5_v2'
    
                    DiffusionBValue = [0, 0, 0, 0,...
                                    3000, 3000, 3000,...
                                    3000, 3000, 3000];
    
                    DiffusionEffBValue = [40.874806178439, 40.874806178439, 40.874806178439, 40.874806178439,...
                                        3004.32428616017, 3001.62509019329, 3001.19552738335,...
                                        3009.2638063984, 3001.13797189685, 2999.56521800681];
    
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    
 
    
                case '40u_DtiStandard_2012'
    
                    DiffusionBValue = [0, 0,...
                                1200, 1200, 1200,...
                                1200, 1200, 1200];
    
                    DiffusionEffBValue = [...
                    40.874806178439, 40.874806178439, 1208.20116109981, 1201.34090122179,...
                    1200.98664722903, 1211.12175725522, 1200.93918208441, 1199.64215559077...    
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0];
    
    
            end
    
        case '20241218_UQ3'

            switch SeriesDescription

                case '40u_verdict_seq1_LR'

                    DiffusionBValue = [0, 0, 600, 600, 600, 600, 600, 600];

                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 615.400011401037, 600.741695507857,... 
                        600.546228372043, 624.156238401682, 603.737798260018, 611.650682773699...  
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0]; 

                case '40u_verdict_seq2_LR'

                    DiffusionBValue = [0, 0, 900, 900, 900, 900, 900, 900];

                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 912.881316083634, 900.806818909437,... 
                        900.594094361001, 923.184038820087, 901.485234531137, 910.699217185179...
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0]; 

                case '40u_verdict_seq3_LR'

                    DiffusionBValue = [0, 0, 1200, 1200, 1200, 1200, 1200, 1200];

                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 1212.58939971525, 1200.85650837969,...
                        1200.63061634199, 1222.36703414999, 1200.60034992734, 1209.90395810655...
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0]; 

                case '40u_verdict_seq4_LR'

                    DiffusionBValue = [0, 0, 1500, 1500, 1500, 1500, 1500, 1500];

                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 1512.42699031062, 1500.88218898693,... 
                        1500.6494917022, 1521.9192709351, 1500.61831347838, 1509.46944383666...
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0]; 

                    

                case '40u_verdict_seq5_LR'

                    DiffusionBValue = [0, 0, 1800, 1800, 1800, 1800, 1800, 1800];

                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 1812.2517664077, 1800.90859884585,... 
                        1800.66890306562, 1821.44065604551, 1800.63678714002, 1809.0058829866...
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0]; 


                case '40u_DtiStandard_2012'

                    DiffusionBValue = [0, 0, 1200, 1200, 1200, 1200, 1200, 1200];

                    DiffusionEffBValue = [...
                        40.874806178439, 40.874806178439, 1212.81394810932, 1200.81874320395,... 
                        1200.60285877035, 1222.99391390378, 1201.04836987655, 1210.51384296977...
                    ];
    
                    % TBC
                    DiffusionDirection = [0, 0, 0, 0, 0, 0, 0, 0]; 
                    
            end


        case '20250131_CPC1'

            switch SeriesDescription

                case '192670001' % STEAM DELTA 80 B 1000

                    DiffusionBValue = [0, 1000, 1000, 1000];
                    
                    % TBC
                    DiffusionEffBValue = [0, 1000, 1000, 1000];
                    DiffusionDirection = [0, 0, 0, 0];

                case '192700001' % STEAM DELTA 120 B 2000

                    DiffusionBValue = [0, 2000, 2000, 2000];

                    % TBC
                    DiffusionEffBValue = [0, 2000, 2000, 2000];
                    DiffusionDirection = [0, 0, 0, 0];  

                case '192720001' % STEAM DELTA 40 B1500

                    DiffusionBValue = [0, 1500, 1500, 1500];

                    % TBC
                    DiffusionEffBValue = [0, 1500, 1500, 1500];
                    DiffusionDirection = [0, 0, 0, 0];    


              case '192690001' % SE DELTA 15 B 2000  

                    DiffusionBValue = [0, 2000, 2000, 2000];

                    % TBC
                    DiffusionEffBValue = [0, 2000, 2000, 2000];
                    DiffusionDirection = [0, 0, 0, 0]; 

              case '192710001' % SE DELTA 20 B 1500  

                    DiffusionBValue = [0, 1500, 1500, 1500];

                    % TBC
                    DiffusionEffBValue = [0, 1500, 1500, 1500];
                    DiffusionDirection = [0, 0, 0, 0];    


              case '192730001' % SE DELTA 20 B 1500  

                    DiffusionBValue = [0, 1000, 1000, 1000];

                    % TBC
                    DiffusionEffBValue = [0, 1000, 1000, 1000];
                    DiffusionDirection = [0, 0, 0, 0];                      
            end



        case '20250224_UQ4'
            
            switch SeriesDescription

                case 'SE_b0_SPOIL10%'

                    DiffusionBValue = [0, 0, 0, 0];
                    DiffusionDirection = [0, 0, 0, 0];
                    DiffusionEffBValue = [40.874806178439 40.874806178439 40.874806178439 40.874806178439];


                % case 'STEAM_DELTA_40'
                % 
                %     DiffusionBValue = [0, 1000, 1000, 1000, 1000, 1000, 1000];
                %     DiffusionDirection = [0,0,0,0,0,0,0];
                %     DiffusionEffBValue = [31.4555326508071 1030.00865891474 1022.785232274 1016.82344802067 1007.11154118745 1016.024651232 994.196831972685];
                % 
                % 
                % case 'STEAM_DELTA_60'
                % 
                %     DiffusionBValue = [0, 1250, 1250, 1250, 1250, 1250, 1250];
                %     DiffusionDirection = [0,0,0,0,0,0,0];
                %     DiffusionEffBValue = [49.3776905745272 1293.67897606535 1283.17531330854 1274.50621838444 1262.69815513832 1273.34467933978 1241.60461141643];
                % 
                % 
                % case 'STEAM_DELTA_80'
                % 
                %     DiffusionBValue = [0, 1500, 1500, 1500, 1500, 1500, 1500];
                %     DiffusionDirection = [0,0,0,0,0,0,0];
                %     DiffusionEffBValue = [67.2998484982473 1556.88949057735 1543.21734229189 1531.93316894712 1518.27178656765 1530.42124549473 1489.10661789428];
                % 
                % case 'STEAM_DELTA_100'
                % 
                %     DiffusionBValue = [0, 1750, 1750, 1750, 1750, 1750, 1750];
                %     DiffusionDirection = [0,0,0,0,0,0,0];
                %     DiffusionEffBValue = [85.2220064219674 1819.8907093611 1803.10094285096 1789.24367391537 1773.83542155106 1787.3869911749 1736.65151519978];
                % 
                % case 'STEAM_DELTA_120'
                % 
                %     DiffusionBValue = [0, 2000, 2000, 2000, 2000, 2000, 2000];
                %     DiffusionDirection = [0,0,0,0,0,0,0];
                %     DiffusionEffBValue = [103.144164345687 2082.77829673617 2062.89852901701 2046.49095792332 2029.39237426483 2044.29256990009 1984.21969891062];
                % 
                % 
                case 'STEAM_DELTA_15'

                    DiffusionBValue = [0, 1000, 1000, 1000, 1000, 1000, 1000];
                    DiffusionDirection = [0,0,0,0,0,0,0];
                    DiffusionEffBValue = [9.05283524615702 1012.54659566948 1009.51440180003 1007.01181036206 1002.41998289654 1006.67649765716 997.513785410845];

                case 'STEAM_DELTA_20'

                    DiffusionBValue = [0, 1250, 1250, 1250, 1250, 1250, 1250];
                    DiffusionDirection = [0,0,0,0,0,0,0];
                    DiffusionEffBValue = [13.533374727087 1269.19965780595 1264.56105207128 1260.73262437321 1253.70811402484 1260.21966791236 1246.20268502774];

                case 'STEAM_DELTA_30'

                    DiffusionBValue = [0, 1500, 1500, 1500, 1500, 1500, 1500];
                    DiffusionDirection = [0,0,0,0,0,0,0];
                    DiffusionEffBValue = [22.4944536889471 1529.76232974936 1522.57766985809 1516.6478813608 1505.76773317781 1515.85337156174 1494.14269803435];    
            
                case 'STEAM_DELTA_40'

                    DiffusionBValue = [0, 1750, 1750, 1750, 1750, 1750, 1750];
                    DiffusionDirection = [0,0,0,0,0,0,0];
                    DiffusionEffBValue = [31.4555326508071 1789.60487926892 1780.04918401636 1772.16248476409 1757.69173988364 1771.10577593863 1742.23028522713];  

                case 'STEAM_DELTA_50'

                    DiffusionBValue = [0, 2000, 2000, 2000, 2000, 2000, 2000];
                    DiffusionDirection = [0,0,0,0,0,0,0];
                    DiffusionEffBValue = [40.4166116126672 2049.11559200975 2037.26951111784 2027.49246456673 2009.55325555165 2026.1824752991 1990.38587549958];  

                case '40u_DtiSTE_SPOIL5%'

                    DiffusionBValue  = [0, 0, 1200, 1200, 1200, 1200, 1200, 1200];
                    DiffusionDirection = [0,0,0,0,0,0,0,0];
                    DiffusionEffBValue = [497.009991353231 497.009991353231 1450.51908514414 1279.40842320093 1232.78246418072 1492.57145401906 1294.52249145437 1339.17154968835];



            end

    end


    % Append to structure
    DTIstruct.DiffusionBValue = DiffusionBValue;
    DTIstruct.DiffusionEffBValue = DiffusionEffBValue;
    DTIstruct.DiffusionDirection = DiffusionDirection;

    % Save structure
    
    folder = fullfile(ImagingDataFolder, ImageType, SampleName, SeriesDescription);
    mkdir(folder);
    save(fullfile(folder, 'DTIstruct.mat'), "DTIstruct");
    
end



