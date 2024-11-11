%{
P1
The Implicit Method
==========================================
Application: PDE of riksy asset (St) of Wang(1996)'s paper
We want to solve: 2nd-order PDE of S with exogenous state var Y
where Y ~ GMB

-I use forward and backward approx: SY
-I change the sign of St (\partial S / \partial t): to start from t0
   : same result: inverse U-shape of S
- Boundaries in S(n+1): the correct to do that is to substitute Eq1 and EqI+1 by
  Boundaries conditions

----------------------------
Author: Hamilton Galindo Gil
Date:   2023 (March, April), 2024 (Nov)
Paper base: Wang(1996) 
----------------------------
%}
%=========================================
clear; clc;
close all;
tic;
%% STEP 1: Parameters
% A. Preferences
% %{
rho    = 0.1;   % impatience rate in discount factor e^(-rho*t)
mu     = 0.05; %E[dY/Y]
sigma  = 0.2;  %Volatility term of dY/Y
lambda = 2/3; %2/3;   % the weight of agent 1(more RRA, RRA=1) in the RA utility function
gamma1 = 1;      % more risk-averse
gamma2 = gamma1/2; % less risk-averse
%}

% longstaff2012  
%rho=0.01; mu=0.03; sigma=0.1; gamma2=2; gamma1=2*gamma2; lambda=1/2;

a = ( (1-lambda)/lambda )^(1/gamma1);
b = 4*a^2; % coming from the agent's 1 budget constraint at t=0                 
                                             
%% STEP 2: Discretization 
% A. State space: structured grid
Ymax = 100;   
Ymin = ( (1 + b/(2+sqrt(b)))^2 - 1 )/b;  %consistent with m0=1
I = 500;                % N of points in the grid: I + 1
deltaY = (Ymax-Ymin)/I; % the distance between grid points
  Y = Ymin:deltaY:Ymax; %the vector of the state variable (grid)
% W = [Wmin... Wmax]

% B. Boundaries of S
r1 = rho + gamma1*mu - (gamma1*(1+gamma1)/2)*sigma^2; % aux variable
r2 = rho + gamma2*mu - (gamma2*(1+gamma2)/2)*sigma^2; % aux variable

N1 = r1 - mu + gamma1*sigma^2;
N2 = r2 - mu + gamma2*sigma^2;

s0_min = Y(1)/N1;    % lower boundary of S
s0_max = Y(end)/N2;  % upper boundary of S

%% STEP 3: Preliminary for iteration of s
Smatrix = [];   %storage of s for every iteration

%% STEP 4: INITIAL GUESS of s (for every point of the state var)
% A. Initial guess of "s"    
    % Sn = [Sn_1, Sn_2, ..., Sn_I]
    s0 = Y.^0.5; %ones(1,length(Y)); %exp(Y); % I tried: exp(Y), it takes 23 iterattion, same results

% B. Initial guess of price function: s0 is a vector
    s = s0;   %S = [S(Y0)  S(Y1) ...        S(YI)]
              %Y   [Y0     Y1        YI-1    YI    ]
        %position  [1      2         I       I+1   ]

%% STEP 5: Iteration of s
maxit= 1000;
crit=10^(-6);  %the criterion to stop iteration and
               %get the solution of "s"
deltat = 1000; %time length (from Achdou et al (2022))

for n=1:maxit
    %% STEP-5.1: Initial point of S function
    S=s;
    Smatrix = [Smatrix; S]; %We save the initial S of every iteration

    %% STEP-5.2: Finite Difference (Forward/Backward Diff Approx & central)
        % A. Forward and Backward Difference
          % S: forward difference (SY)                  
            dSf  = [(S(2:end) - S(1:end-1))/deltaY 0];
            % Boundary nodes (Ymax): dSf(I+1,:) = 0
                             % we will use it since we do not backward
                             % ghost node: S(I+2)
          
          % S: backward difference (SY)                        
            dSb  = [0 (S(2:end) - S(1:end-1))/deltaY];
            % Boundary nodes (Ymin): dSb(1,:) = 0
                             % it will never be used
                             % because at Ymin we use forward
                             % ghost node: S(0)

          % Central difference (SYY)
             ddSYY = [(S(2) - S(1))/deltaY^2,...
                    (S(3:end) - 2*S(2:end-1) + S(1:end-2))/deltaY^2,...
                    (-S(end) + S(end-1))/deltaY^2 ];
                          
    %% STEP-5.3: Upwind scheme                    
       % A. Implementation of "eta", "psi", and "r"
         %aux
         eta = 2*(1+b*Y - sqrt(1+b*Y));
         %psi (- Sharpe ratio)
         psi = -sigma*gamma1*b*Y./eta;
         %r
         r = rho + (mu.*Y).*(gamma1*b./eta)...
             - (1/2).*( (sigma*Y).^2 ).*(gamma1*b).*( ((2+gamma1)*b - b.*(1+b.*Y).^(-0.5) )./(eta.^2) );

       % B. Implementation of Upwind Scheme 
         % (B.1) bi: coefficient
            bcoef = Y.*(mu + psi.*sigma); %row vector

         % (B.2) Indicator Functions
            % dS_upwind makes a choice of forward or backward differences based on
            % based on the sign of the drift (Sy):
            If = bcoef > 0; %positive drift --> forward difference
            Ib = bcoef < 0; %negative drift --> backward difference: 
                          %Ib is a logic vector: zeros and ones: 
                          %1 means "true"
            I0 = (1-If-Ib); %when b=0 

         % (B.3) Boundaries conditions   
          % To be sure that in Ymin we will use Forward
            % consistent with X1=0 
             %If(1)=1; Ib(1)=0; I0(1)=0;
          % To be sure that in Y(I+1) we will use Backward
            % consistent with Z(I+1)=0
             %If(end)=0; Ib(end)=1; I0(end)=0;
          % Already taken care of automatically
                  
         % (B.4) The first derivative with Upwind scheme            
            SY_Upwind = dSf.*If + dSb.*Ib + dSf.*I0;
            
            %storage of "b" for every "i (grid)"
                b_Upwind  = bcoef.*If  + bcoef.*Ib;      
      
            %check
                check = [Y'     dSf'    dSb' ...
                         ddSYY' bcoef'  If'...
                         Ib'    I0'     SY_Upwind'];         
         
    %% STEP-5.4: Discretization of Price Equation
       % Implicit method: 
       % We need to construct a matrix

       % A. Coefficients (column vectors)        
        X = - min(bcoef',0)/deltaY + (1/2)*( (sigma.*Y').^2 )/deltaY^2;
        H =   min(bcoef',0)/deltaY - max(bcoef',0)/deltaY - ((sigma.*Y').^2)/deltaY^2; 
        Z =   max(bcoef',0)/deltaY + (1/2)*((sigma.*Y').^2)/deltaY^2;
        
       % B. Matrix of coefficients: "A"        
        % Up Diagonal (Z)
            updiag = [ 0; 0; Z(2:end-1)]; % spdiags counts since the 2nd position
        % Central Diagonal (H)
            centdiag = [-1; H(2:end-1); -1];
        % Down Diagonal (X)
            lowdiag = [ X(2:end-1); 0; 0 ];    
        % An    
            An = spdiags([lowdiag centdiag updiag], -1:1, I+1, I+1);
            
        % See the diagonal matrix XYZ
            %spy(An)
        
            %{ 
            look at this example to undertand how "spdiags" works
            ZZL = [1 2 3]'
            ZZC = [-1 -1 -1]'
            ZZU = [4 5 6]'
            
            ZZZ = spdiags([ZZL ZZC ZZU], -1:1, 3,3)
            %}

       % C. Left-hand matrix: "B" 
        R = diag([0 r(2:end-1) 0]);
        bc = [0; 1/deltat*ones(size(Y',1)-2,1); 0]; %boundary conditions  
        B = diag(bc) + R - An;

       % D. Right-hand matrix: "bn" (column vector)       
        Ytilde = [s0_min Y(2:end-1) s0_max]'; 
        bn = Ytilde + diag(bc)*S';
        
       % E. Solve the system of equations: finding S^n+1        
        S = B\bn; % column vector 

            %S(1)   = Y(1)/rho;
            %S(end) = Y(end)/bS_upper; 

    %% STEP-5.5: Update of the value function        
        Schange = S - s';  % since we have "S", we calculate "Schange"
        s = S';            % the "new initial guess"

    %% STEP-5.6: The optimal value function
        %We use the "sup-norm" 
        %We can use others: e.g., Euclidean norm        
        dist(n) = max(abs(Schange));

        %% Stability
        imax = I+1;
        error(n) = (1/(imax*deltat))*sqrt(sum(Schange.^2));

        if dist(n)<crit %crit=10^(-6)
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end
        
    % To know in what "iteration" we are
    disp(n)

end
toc;

%% Policy functions
t=0;      %Analysis at t=0
T = 100;

% Optimal consumption (Pareto optimality) and Sharing rule
c1 = (2/b)*( sqrt(1 + b*Y) - 1 );
c2 = Y - c1;

% The Stochastic Discount Factor
% Eq: dm = -r*m*dt + psi*m*dZ
m = ( (b/2)^(gamma1)/( 1 + (sqrt(b)/2)^(gamma1) ) )*...
    ( (exp(-rho*t))./(sqrt(1+b.*Y) - 1).^(gamma1) );

% Wealth (t=0, T=100)
    % See Lemma 3.10
betas = [];
gs = [];
gsYt = [];
A = ( (b/2)^(gamma1-1) )/(1 + (sqrt(b)/2)^(gamma1));
for i = 1:length(Y)
Yt = Y(i);
n1 = 401;
h = (T - t)/(n1-1);
    for ii = 1:n1
        s = t + (ii-1)*h;
        st = s-t;
        gs(ii)    = h_fun(st,b,mu,sigma,gamma1,Yt);
        betas(ii) = exp(-rho*s);
    end
omega_tilde = [h/2,h*ones(1,n1-2) ,h/2];    
gsYt(i) = sum(omega_tilde.*betas.*gs); 
end
W1 = (A./m).*gsYt; % wealth of agent 1
W2 = S' - W1;      % wealth of agent 2

% ----Special case: gamma1=1, gamma2=0.5 (analytic)
% Wealth (t=0, T=100)
% See Lemma 3.10 
gt = -( 2/(rho*(2+sqrt(b))) )*(exp(-rho*T) - exp(-rho*t));
W1old = gt.*m.^(-1); % wealth of agent 1
W2old = S' - W1;      % wealth of agent 2

    % Graph
    figure('Name','W1')
    subplot(2,2,1)
    plot(Y,W1,'--',Y,W1old,':','LineWidth', 1.5)
    legend('Numerical approximation','Analytical solution')
    title('Wealth of the More Risk-Averse Agent ($W_{1t}$)')
    xlabel('Endowment (Y)')
    grid;
    % Remark: W1old = W1 : our code to calculate expectations is good!
    % Save the figure
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);
    print(h, '-dpdf', strcat('W1t','_Fig0.pdf'));


% Stock volatility
dSY = [0 (S(2:end) - S(1:end-1))'/deltaY]; % 1st derivative of S
Compare = [SY_Upwind' dSY']; % they are almost the same except in I+1
sigmat = sigma*Y.*dSY./S';
%sigmat = sigma*Y.*SY_Upwind./S';
plot(dSY)
title('dSY (slope of S)')

% Optimal portfolio: agent 1
% See Lemma 3.10 
    % possition=2 (backward diff)
    fy = ( gsYt(2:end) - gsYt(1:end-1) )/deltaY;  % df/dY
    my = ( m(2:end) - m(1:end-1) )/deltaY;

    % All positions
    fy1 = [fy 0];
    my1 = [my 0];

% Optimal portfolio: agent 1 (more risk-averse)   
w11 = sigma*Y.*(fy1./gsYt - my1./m)./sigmat; % risky asset
w12 = 1 - w11;                      % riskless asset

% N of shares
N11 = w11.*W1./S'; % fraction of 1share held by agent1
N21 = 1 - N11;     % fraction of 1share held by agent2

% Optimal portfolio: agent 2 (less risk-averse)
w21 = S'.*N21./W2;   % risky asset
w22 = 1 - w21;       % riskless asset

% NB = amount of money invested in the riskless asset
NB1 = W1 - N11.*S'; % agent1 (more risk-averse)
                    % NB1>0 : lender : buy riskless assets
NB2 = W2 - N21.*S'; % agent2 (less risk-averse)
                    % NB2<0 : borrower : sell riskless assets (leverage)

% Expected rate of return
beta = r - psi.*sigmat;

% Price-Dividend Ratio
pd = S./Y';
w  = c1./Y;   % relative consumption of agent 1 (more risk-averse)
s  = c2./Y;   % relative consumption of agent 2 (less risk-averse): longstaff2012

% choose the beginning and end of the graph
tend = 500;
tinit = 2;

