clear;

delta = 2;
DELTAS = [15, 20, 30, 40, 50];
% DELTAS = [40,60,80,100,120];
bvals = [1000,1250,1500,1750,2000];

Rs = [5,10,20,30,40,50];
d = 1;

figure
for R = Rs

    SPHEREsignals = zeros(size(DELTAS));
    BALLsignals = zeros(size(DELTAS));

    for Dindx = 1:length(DELTAS)
    
        DELTA = DELTAS(Dindx);
        bval = bvals(Dindx);
        G = stejskal(delta, DELTA, bval=bval);
        SPHEREsignals(Dindx) = sphereGPD(delta, DELTA, G, R, d);
        BALLsignals(Dindx) = ball(bval, d);

    end

    plot(bvals, SPHEREsignals, '--', DisplayName=['R=' num2str(R)])
    hold on
end
plot(bvals, BALLsignals, '-', color='k', DisplayName='ball')
legend
ylim([0,1])