%{
The Implicit Method
==========================================
Application: PDE of riksy asset (St) of Wang(1996)'s paper
We want to solve: 2nd-order PDE of S with exogenous state var Y

I use forward and backward approx: SY
- S has an inverse U-shape

- This is the same m-file: ChXX_PDE_S2A.m
----------------------------
Author: Hamilton Galindo Gil
Date:   2023 (March)
Paper base: Wang(1996) 
----------------------------
Book: Heterogeneous Agents in Asset Pricing
Chapter: XX
%}
%=========================================
clear; clc;
close all;
tic;
%% STEP 1: Parameters
% A. Preferences
rho    = 0.1;   % impatience rate in discount factor e^(-rho*t)
lambda = 2/3;   % the weight of agent 1(more RRA, RRA=1) in the RA utility function
b = 4*( (1-lambda)^2 )/lambda^2; % coming from the agent's 1 budget constraint at t=0                 
                 
% B. Exogeneous State Variable Dynamic (Y)
mu      = 0.05; %E[dY/Y]
sigma   = 0.3;  %Volatility term of dY/Y
                                             
%% STEP 2: Discretization 
% A. State space: structured grid
Ymax = 100;   
Ymin = ( (1 + b/(2+sqrt(b)))^2 - 1 )/b;  %consistent with m0=1
I = 500;                % N of points in the grid: I + 1
deltaY = (Ymax-Ymin)/I; % the distance between grid points
  Y = Ymin:deltaY:Ymax; %the vector of the state variable (grid)
% W = [Wmin... Wmax]

%% STEP 3: Preliminary for iteration of s
Smatrix = [];   %storage of s for every iteration

%% STEP 4: INITIAL GUESS of s (for every point of the state var)
% A. Initial guess of "s"
    % Sn = [Sn_1, Sn_2, ..., Sn_I]
    s0 = Y.^0.5; % I tried: exp(Y), it takes 23 iterattion, same results

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
    %% STEP-5.1: Initial point of value function
    S=s;
    Smatrix = [Smatrix; S]; %We save the initial S of every iteration

    %% STEP-5.2: Finite Difference (Forward/Backward Diff Approx & central)
        % A. Forward and Backward Difference
          % Boundaries
                r2 = rho + 0.5*mu - ((1+0.5)/2)*sigma^2;
                bS_upper = r2 - mu + 0.5*sigma^2;

            Slower = (Ymin-deltaY)/rho;
            Supper = (Ymax+deltaY)/bS_upper;

          % S: forward difference (SY)                  
            %dSf  = [(S(2:end) - S(1:end-1))/deltaY 0];
            dSf  = [(S(2:end) - S(1:end-1))/deltaY (Supper - S(end))/deltaY];
            % Boundary nodes (Ymax): dSf(I+1,:) = 0
                             % we will use it since we do not backward
                             % ghost node: S(I+2)
          
          % S: backward difference (SY)                        
            %dSb  = [0 (S(2:end) - S(1:end-1))/deltaY];
            dSb  = [(S(1) - Slower)/deltaY (S(2:end) - S(1:end-1))/deltaY];
            % Boundary nodes (Ymin): dSb(1,:) = 0
                             % it will never be used
                             % because at Ymin we use forward
                             % ghost node: S(0)

          % Central difference (SYY)
             %ddSYY = [(S(2) - S(1))/deltaY^2,...
             %       (S(3:end) - 2*S(2:end-1) + S(1:end-2))/deltaY^2,...
             %       (-S(end) + S(end-1))/deltaY^2 ];
             ddSYY = [(S(2) - 2*S(1) + Slower)/deltaY^2,...
                    (S(3:end) - 2*S(2:end-1) + S(1:end-2))/deltaY^2,...
                    (Supper -2*S(end) + S(end-1))/deltaY^2 ];
                          
    %% STEP-5.3: Upwind scheme                    
       % A. Implementation of "eta", "psi", and "r"
         %aux
         eta = 2*(1+b*Y - sqrt(1+b*Y));
         %psi (- Sharpe ratio)
         psi = -sigma*b*Y./eta;
         %r
         r = rho + (mu.*Y).*(b./eta)...
             - (1/2).*( (sigma*Y).^2 ).*(b).*( (3*b - b.*(1+b.*Y).^(-0.5) )./(eta.^2) );

       % B. Implementation of Upwind Scheme 
         % (B.1) bi: coefficient
            bcoef = -Y.*(mu + psi.*sigma); %row vector

         % (B.2) Indicator Functions
            % dS_upwind makes a choice of forward or backward differences based on
            % based on the sign of the drift (capital):
            If = bcoef > 0; %positive drift --> forward difference
            Ib = bcoef < 0; %negative drift --> backward difference: 
                          %Ib is a logic vector: zeros and ones: 
                          %1 means "true"
            I0 = (1-If-Ib); %when b=0 

         % (B.3) Boundaries conditions   
          % To be sure that in Wmin we will use Backward
             If(1)=0; Ib(1)=1; I0(1)=0;
          % To be sure that in W(I+1) we will use Forward
             If(end)=1; Ib(end)=0; I0(end)=0;
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
        X = - min(bcoef',0)/deltaY -(1/2)*( (sigma.*Y').^2 )/deltaY^2;
        H =   min(bcoef',0)/deltaY - max(bcoef',0)/deltaY + ((sigma.*Y').^2)/deltaY^2; 
        Z =   max(bcoef',0)/deltaY - (1/2)*((sigma.*Y').^2)/deltaY^2;
        
       % B. Matrix of coefficients: "A"        
        % Up Diagonal (Z)
            updiag = [ 0; Z(1:end-1)]; % spdiags counts since the 2nd position
        % Central Diagonal (H)
            centdiag = H;
        % Down Diagonal (X)
            lowdiag = [ X(2:end); 0 ];    
        % An    
            An = spdiags([lowdiag centdiag updiag], -1:1, I+1, I+1);
            
        % See the diagonal matrix XYZ
            spy(An)
        
            %{ 
            look at this example to undertand how "spdiags" works
            ZZL = [1 2 3]'
            ZZC = [-1 -1 -1]'
            ZZU = [4 5 6]'
            
            ZZZ = spdiags([ZZL ZZC ZZU], -1:1, 3,3)
            %}

       % C. Left hand matrix: "B" 
        R = diag(r);
        B = (1/deltat)*speye(I+1) - R - An;
        
       % D. Right-hand matrix: "bn" (column vector)
        bn = -Y' + (1/deltat)*S' ;

       % E. Solve the system of equation: finding S^n+1        
        S = B\bn; % column vector 
            
    %% STEP-5.5: Update of the value function        
        Schange = S - s';  % since we have "S", we calculate "Schange"
        s = S';            % the "new initial guess"

    %% STEP-5.6: The optimal value function
        %We use the "Absolute-value norm" 
        %We can use others: e.g., Euclidean norm        
        dist(n) = max(abs(Schange));
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
m = (b/(2 + sqrt(b)))*( (exp(-rho*t))./(sqrt(1+b.*Y) - 1) );

% Stock volatility
sigmat = sigma*Y.*dSf./S';

% Optimal portfolio: agent 1
w11 = -psi./sigmat; % risky asset in his portfolio

% Wealth (t=0, T=100) 
gt = -( 2/(rho*(2+sqrt(b))) )*(exp(-rho*T) - exp(-rho*t));
W1 = gt.*m.^(-1); % wealth of agent 1
W2 = S' - W1;      % wealth of agent 2

% N of shares
N11 = w11.*W1./S'; % fraction of 1share held by agent1
N21 = 1 - N11; % fraction of 1share held by agent2

% Optimal portfolio: agent 2
w21 = S'.*N21./W2;

%NB = amount of money invested in the riskless asset
NB1 = W1 - N11.*S'; %agent1 (RRA = 1)
                    % NB1>0 : lender : buy riskless assets
NB2 = W2 - N21.*S'; %agent2 (RRA = 1/2)
                    % NB2<0 : borrower : sell riskless assets (leverage)
% expected rate of return
beta = r - psi.*sigmat;

% Price-Dividend Ratio
pd = S./Y';
w = c1./Y;   % relative consumption

%% Graphs (policy function)

% Plot (Fig 0)
figure('Name','Asset Price: S')
yyaxis left
plot(Y,S)

yyaxis right
plot(Y, bcoef)

legend('S', 'bcoef')
title('S')
xlabel('Endowment (Y)')


% choose the end
tend = 200;
tinit = 10;

% Plot (Fig 1)
figure('Name','Fig1')
subplot(2,2,1)
plot(Y(tinit:tend),c1(tinit:tend),'r:',...
     Y(tinit:tend),c2(tinit:tend),'b--','LineWidth', 1.5);
xlabel('Endowment (Y)')
titlestr = strcat('Endowment and Optimal Consumption','($\lambda$=',num2str(lambda),')');
title(titlestr,'interpreter','latex')
legend('c_1 (RRA=1)', 'c_2 (RRA=1/2)')
grid;

subplot(2,2,2)
plot(Y(tinit:tend),S(tinit:tend),'k',...
     Y(tinit:tend),W1(tinit:tend),'r:',...
     Y(tinit:tend),W2(tinit:tend),'b--','LineWidth', 1.5);
xlabel('Endowment (Y)')
title('Asset Price and Wealth')
legend('Asset Price (S)','Wealth of Agent 1','Wealth of Agent 2')
grid;

subplot(2,2,3)
plot(Y(tinit:tend),r(tinit:tend),'k',...
     Y(tinit:tend),-psi(tinit:tend),'r:',...
     Y(tinit:tend),beta(tinit:tend),'b--','LineWidth', 1.5);
xlabel('Endowment (Y)')
title('Asset Prices I')
legend('Interest Rate (r)', 'Price of Risk (-\psi)', 'Expected Rate of Return (\beta)')
grid;

subplot(2,2,4)
plot(Y(tinit:tend),m(tinit:tend),'r:',...
     Y(tinit:tend),sigmat(tinit:tend),'b--','LineWidth', 1.5);
xlabel('Endowment (Y)')
title('Asset Prices II')
legend('Stochastic Discount Factor (m)', 'Stock Volatility (\sigma_t)')
grid;


% Plot (Fig 2: N of shares)
figure('Name','Fig2')
subplot(2,2,1)
plot(Y(tinit:tend),N11(tinit:tend),'r:',...
     Y(tinit:tend),N21(tinit:tend),'b--',...
    'LineWidth', 1.5);
xlabel('Endowment (Y)')
titlestr1 = strcat('Riksy Asset Shares', '($\lambda$=',num2str(lambda),')');
title(titlestr1,'Interpreter','latex')
legend('N_1^{(1)}', 'N_2^{(1)}')
grid;

subplot(2,2,2)
plot(Y(tinit:tend),NB1(tinit:tend),'r:',...
     Y(tinit:tend),NB2(tinit:tend),'b--',...
    'LineWidth', 1.5);
xlabel('Endowment (Y)')
title('Money invested in riskless asset')
legend('B*N_1^{(2)}', 'B*N_2^{(2)}')
grid;

subplot(2,2,3)
plot(Y(tinit:tend), w11(tinit:tend),'r:',......
     Y(tinit:tend), w21(tinit:tend),'b--','LineWidth', 1.5);
xlabel('Endowment (Y)')
title('Optimal Portfolio: risky asset')
legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;

subplot(2,2,4)
%plot(Y(tinit:tend), pd(tinit:tend),'LineWidth', 1.5);
plot(Y, pd,'LineWidth', 1.5);
xlabel('Endowment (Y)')
title('Price-Dividend Ratio (S/Y)')
%legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;

%Figure
figure('Name','Fig3')
%subplot(2,2,1)
plot(w(tinit:tend), r(tinit:tend),'LineWidth', 1.5);
xlabel('Relative Consumption of Agent 1 (w = c_1/Y)')
title('r')
%legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;


%% Saving "S"
Wrong_S = S;
save('Wrong_S.mat','Wrong_S')   % save variable in the output.mat file