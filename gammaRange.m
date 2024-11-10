function results = gammaRange(rho,mu,sigma)
%{ 
This rule for gamma1 comes from:
    (1) rk > 0
    (2) Nk > 0
%}
aux1 = sqrt( (0.5*sigma^2 - mu)^2 + 2*rho*sigma^2 ); 
R1r = (1/sigma^2)*( mu - 0.5*sigma^2 - aux1 );
R2r = (1/sigma^2)*( mu - 0.5*sigma^2 + aux1 );

aux2 = sqrt( (0.5*sigma^2 + mu)^2 - 2*(mu - rho)*sigma^2 ); 
R1N = (1/sigma^2)*( mu + 0.5*sigma^2 - aux2 ); 
R2N = (1/sigma^2)*( mu + 0.5*sigma^2 + aux2 );

gammakMin = max(R1r,R1N);
gammakMax = min(R2r,R2N);

gRange = [gammakMin,gammakMax];
Max_gamma1 = 2*(2*mu/sigma^2 + 1)/3;   % Max gamma1

results.Rr = [R1r, R2r];
results.RN = [R1N, R2N];
results.gRange = gRange;
results.Max_gamma1 = Max_gamma1;

end