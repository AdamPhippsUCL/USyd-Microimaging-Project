% MFP imaging idea

Rs = linspace(5,10,11);
D = 2;



deltas = linspace(5, 25, 21);
DELTAs = linspace(15, 50, 36);
bval = 1000;
Gmax = 65;

signals = zeros(length(Rs), length(deltas), length(DELTAs));

for Rindx = 1:length(Rs)
    for dindx = 1:length(deltas)
        for Dindx = 1:length(DELTAs)
    
            R = Rs(Rindx);
            delta = deltas(dindx);
            DELTA = DELTAs(Dindx);
            G = stejskal(delta, DELTA, bval = bval);

            if G<Gmax
                signals(Rindx, dindx, Dindx) = sphereGPD(delta, DELTA, G, R, D);
            end
        
        end
    end
end

% figure;
% tiledlayout(1,1);
% nexttile;
% imshow(squeeze(signals(:,1,:)), [])
% axis('on', 'image');
% xticklabels(DELTAs(xticks))
% yticklabels(Rs(yticks))



signalrange = squeeze(range(signals, 1));
figure;
tiledlayout(1,1);
nexttile;
imshow(signalrange, [min(signalrange(signalrange>0)) max(signalrange(signalrange>0))])
colorbar;
axis('on', 'image');
xticklabels(DELTAs(xticks))
yticklabels(deltas(yticks))







% scheme = load("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\VERDICT-Processing\Schemes\Original.mat").scheme;
% 
% 
% figure;
% 
% 
% for Rindx = 1:length(Rs)
% 
%     R = Rs(Rindx);
% 
%     signals = zeros(1,length(scheme));
% 
%     for ischeme = 1:length(scheme)
% 
%         delta = scheme(ischeme).delta;
%         DELTA = scheme(ischeme).DELTA;
%         b = scheme(ischeme).bval;
%         G = scheme(ischeme).G;
% 
%         signals(ischeme) = sphereGPD(delta, DELTA, G, R, D);
% 
%     end
% 
%     [sortedb,I] = sort([scheme(:).bval]);
%     plot(sortedb, signals(I));
%     hold on
% 
% 
% 
% end