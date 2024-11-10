function [g] = h_fun(st,b,mu,sigma,gamma1,Yt)

mu_s    = (mu - 0.5*sigma^2)*st; 
sigma_s = sigma*sqrt(st);

% A. Function to be integrated
    f = @(x) ( sqrt(1 + b*Yt.*exp(mu_s + sigma_s.*x)) - 1 ).^(1-gamma1);
    
% D. Distribution parameters
    % x ~ N(0,1)
    mu_x    = 0;          % mean vector
    sigma_x = 1;          % variance matrix

% E. Gaussian integration
    nn = 21;                     % order of approximation
    
    [x,w] = qnwnorm(nn,mu_x,sigma_x);     % Gaussian normal nodes and weights
    fexp = w'*f(x);                  % Gaussian integration of f
    fprintf('Guassian Quadrature:      %10.3f\n',fexp)

% F. g
    g = fexp;
    
