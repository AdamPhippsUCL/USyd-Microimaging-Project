% MATLAB Script to simulate VERDICT contrast for different sequence parameters 

% IC contrast = IC signal - EES signal
% R contrast = std(IC signal) (across range of R values)

% Total contrast = IC contrast * R contrast  (Can change this to give more
% weighting to either R or IC contrast...)

%% Settings

% Noise settings
SNR = 50; % SNR at TE=0;
sigma0 = 1/SNR;

% Tissue settings
dIC=2;
dEES=2;
Rs = linspace(6,9,3); % Different cell sizes to consider
T2 = 30;

% Imaging settings
bvals = [800,1000,1200];%linspace(400,2000,17);
Gmax = 3000;
Deltas = linspace(1,50,51);
deltas = linspace(1,30,31);
TEconst = 5; % TE = delta + DELTA + TEconst



% == Initialise empty arrays

% Maximum IC contrast
maxICcontrast = zeros(size(bvals));
ICbestdeltas = zeros(size(bvals));
ICbestDELTAs = zeros(size(bvals));

% Maximum R contrast
maxRcontrast = zeros(size(bvals));
Rbestdeltas = zeros(size(bvals));
RbestDELTAs = zeros(size(bvals));

% Maximum 'total' contrast
maxtotalcontrast = zeros(size(bvals));
totalbestdeltas = zeros(size(bvals));
totalbestDELTAs = zeros(size(bvals));


%% Simulation

for bindx = 1:length(bvals)

    bval = bvals(bindx);
    disp(['b value: ' num2str(bval) ])
    
    % Initialise arrays for results
    Gs =  zeros(length(deltas), length(Deltas));
    EESSignals = zeros(length(Rs),length(deltas), length(Deltas));
    ICSignals = zeros(length(Rs),length(deltas), length(Deltas));
    noise = zeros(length(Rs),length(deltas), length(Deltas));
    
    for Rindx = 1:length(Rs)
        R = Rs(Rindx);
        for Dindx = 1:length(Deltas)
            for dindx = 1:length(deltas)
                
                % Gradient timings
                delta = deltas(dindx);
                Delta = Deltas(Dindx);
    
                % Echo time
                TE = delta+Delta+TEconst;
    
                % b=0 signal
                b0 = exp(-TE/T2);

                noise(Rindx, dindx, Dindx) = sigma0/b0;
    
                % IC signal
                if delta<Delta
                    G=stejskal(delta, Delta, bval=bval);
                    Gs(dindx, Dindx) = G;
                    spheresignal = sphereGPD(delta,Delta,G,R, dIC);
                    spheredist = makedist('Rician', s=spheresignal, sigma = sigma0*sqrt(2) ); % sqrt(2) as simple way to model normalization by noisy b0 signal
                    ICSignals(Rindx, dindx, Dindx) = (spheredist.mean);

                    % [spheredist, signals] = RatioDist(b0, b0*spheresignal, sigma0, N0=N0, Nb=Nb);
                    % ICSignals(Rindx, dindx, Dindx)=mean(spheredist.*transpose(signals));

                
                end
        
                % EES signal
                normsignal = ball(bval, dEES );
                normdist = makedist('Rician', s=normsignal, sigma = sigma0*sqrt(2));
                EESSignals(Rindx, dindx, Dindx) = (normdist.mean)*b0;

                % [normdist, signals] = RatioDist(b0, b0*normsignal, sigma0, N0=N0, Nb=Nb);
                % EESSignals(Rindx, dindx, Dindx)=mean(normdist.*transpose(signals));
                % 
            
            end
        end
    end
    
    % Remove infinities
    ICSignals(isinf(ICSignals)) = 0;
    ICSignals(ICSignals>1) = 0;
    

    % == Results
    
    % Contrast (difference between IC and EES signals)
    Contrast = (ICSignals-EESSignals)./noise;
    ICcontrast = squeeze(mean(Contrast,1)).*(Gs<Gmax);
    Rcontrast = squeeze(std(Contrast,1)).*(Gs<Gmax);

    % total
    totalcontrast = ICcontrast.*Rcontrast;
    % 
    % 
    f=figure;
    tiledlayout(1,3);
    sgtitle(['b value: ' num2str(bval)])

    % IC contrast
    nexttile;
    imshow(ICcontrast,[min(ICcontrast(ICcontrast>0.001)) max(ICcontrast(:))])
    colorbar;
    xlabel('Delta')
    ylabel('delta')
    title(['IC contrast to noise ratio, Max = ' num2str(max(ICcontrast(:))) ])
    axis('on', 'image');
    xticklabels(Deltas(xticks))
    yticklabels(deltas(yticks))

    % R contrast
    nexttile;
    imshow(Rcontrast,[min(Rcontrast(Rcontrast>0.001)) max(Rcontrast(:))])
    colorbar;
    xlabel('Delta')
    ylabel('delta')
    title(['R contrast to noise ratio, Max = ' num2str(max(Rcontrast(:))) ])
    axis('on', 'image');
    xticklabels(Deltas(xticks))
    yticklabels(deltas(yticks))

    % R contrast
    nexttile;
    imshow(totalcontrast,[min(totalcontrast(totalcontrast>0.001)) max(totalcontrast(:))])
    colorbar;
    xlabel('Delta')
    ylabel('delta')
    title(['Total contrast to noise ratio, Max = ' num2str(max(totalcontrast(:))) ])
    axis('on', 'image');
    xticklabels(Deltas(xticks))
    yticklabels(deltas(yticks))


    % Maximum IC contrast
    [maxICcontrast(bindx),I] = max(ICcontrast(:));
    [I1, I2] = ind2sub(size(ICcontrast),I);
    ICbestdeltas(bindx) = deltas(I1);
    ICbestDELTAs(bindx) = Deltas(I2);


    % Maximum R contrast
    [maxRcontrast(bindx),I] = max(Rcontrast(:));
    [I1, I2] = ind2sub(size(Rcontrast),I);
    Rbestdeltas(bindx) = deltas(I1);
    RbestDELTAs(bindx) = Deltas(I2);


    % Maximum 'total' contrast
    [maxtotalcontrast(bindx),I] = max(totalcontrast(:));
    [I1, I2] = ind2sub(size(totalcontrast),I);
    totalbestdeltas(bindx) = deltas(I1);
    totalbestDELTAs(bindx) = Deltas(I2);

end

% Maximum IC contrast
figure;
tiledlayout(1,2);
sgtitle('IC contrast')

nexttile;
plot(bvals, maxICcontrast, '*-');
xlabel('b value (s/mm^2)')
ylabel('Max IC contrast')

nexttile;
plot(bvals, ICbestdeltas, '*-', DisplayName = 'delta');
hold on
plot(bvals, ICbestDELTAs, '*-', DisplayName = 'DELTA');
legend;
xlabel('b value (s/mm^2)')
ylabel('Best delta/DELTA (ms)')


% Maximum R contrast
figure;
tiledlayout(1,2);
sgtitle('R contrast')

nexttile;
plot(bvals, maxRcontrast, '*-');
xlabel('b value (s/mm^2)')
ylabel('Max R contrast')

nexttile;
plot(bvals, Rbestdeltas, '*-', DisplayName = 'delta');
hold on
plot(bvals, RbestDELTAs, '*-', DisplayName = 'DELTA');
legend;
xlabel('b value (s/mm^2)')
ylabel('Best delta/DELTA (ms)')


% Maximum total contrast
figure;
tiledlayout(1,2);
sgtitle('Total contrast')

nexttile;
plot(bvals, maxtotalcontrast, '*-');
xlabel('b value (s/mm^2)')
ylabel('Max total contrast')

nexttile;
plot(bvals, totalbestdeltas, '*-', DisplayName = 'delta');
hold on
plot(bvals, totalbestDELTAs, '*-', DisplayName = 'DELTA');
legend;
xlabel('b value (s/mm^2)')
ylabel('Best delta/DELTA (ms)')