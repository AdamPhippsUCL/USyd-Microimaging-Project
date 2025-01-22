% Script to save DTI sequence information as structure

%% Initial details

% Imaging data folder
ImagingDataFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\USyd Microimaging Project\USyd-Microimaging-Project\Imaging Data";

% Image type
ImageType = 'MAT';

% Sample name
SampleName = '20241218_UQ3';

% Series Description
SeriesDescriptions = {
    '40u_verdict_seq1_LR',...
    '40u_verdict_seq2_LR',...    
    '40u_verdict_seq3_LR',...
    '40u_verdict_seq4_LR',...
    '40u_verdict_seq5_LR',...
    '40u_DtiStandard_2012',...
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



