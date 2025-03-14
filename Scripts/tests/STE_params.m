
% b
bval = 1500;

% Delta
DELTA = 55;

% delta
delta = 11;

% Gmax
Gmax = 65;

% T2
T2 = 75;


%% SE signal

TE = DELTA + delta + 22;
SE = exp(-TE/T2)

%% STE signal
TM = 40;
TE = DELTA + delta + 22-TM;

STE = (0.5)*exp(-TE/T2)

