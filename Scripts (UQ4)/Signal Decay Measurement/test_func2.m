function signals = test_func2(b,X)

% b = [Ss, Sg, Sl]
% X = [fs1, fg1, fl1; fs2, fg2, fl2; ... ]

Ss = b(1);
Sg = b(2);
Sl = b(3);

N = size(X,1);
signals = zeros(N,1);

for indx = 1:N

    fs = X(indx,1);
    fg = X(indx,2);
    fl = X(indx,3);

    signals(indx) = fs*Ss + fg*Sg + fl*Sl;
    
end


end