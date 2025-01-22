% Function to fit Gaussian distribution(s) to histogram

function coeffs = fitGaussians(counts, bincentres, opts)

arguments
    counts % Histogram bin counts
    bincentres % Histogram bin centres
    
    % Options
    opts.N = 1 % Number of Gaussian distributions
    opts.beta0guess = [1,1,1] % first guesses [weight, mu, sigma] for each distribution
    opts.lb = [0, 0, 0] % Lower bounds for distribution parameters
    opts.ub = [1,2,2] % Upper bounds for distribution parameters


end


% Normalise histogram counts
freqs = counts/sum(counts(:));


% Save N to base
assignin('base', 'thisN', opts.N);


[coeffs, resnorm] = lsqcurvefit(@evalbinfreqs, opts.beta0guess, bincentres, freqs, opts.lb, opts.ub);



end




function binfreqs = evalbinfreqs(b, x)

arguments
    b % [w1, mu1, sigma1, f2, mu2, sigma2, ... ]
    x % bincentres
end

% Number of Gaussians
N = evalin('base', 'thisN');



% Generate pdf for each Gaussian
pdfs = zeros(N, length(x));

for n = 1:N
    w = b(3*n-2);
    mu = b(3*n-1);
    sigma = b(3*n);
    pdfs(n,:) =  w*normpdf(x, mu, sigma);
end

% Normalize
binfreqs = sum(pdfs/sum(pdfs(:)), 1);

end