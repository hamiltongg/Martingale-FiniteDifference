%{
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
Date:   2023 (March, April)
Paper base: Wang(1996) 
----------------------------
Book: Heterogeneous Agents in Asset Pricing
Chapter: XX
%}
%=========================================
clear; clc;
close all;

%% Run the main m-file

run G_Main_PDE_S2F.m

%% Graphs 1 (Solution: asset price)

% Plot (Fig 0)
figure('Name','Asset Price: S')

subplot(2,2,1)
    plot(Y,S,'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    subtitle('(considering boundaries)')
    xlabel('Endowment (Y)')
    grid;

% Load "S" w/o being careful the boundaries of S    
load('Wrong_S.mat'); % This var comes from "Wrong_PDE_S2A.m"

subplot(2,2,2)
    plot(Y,Wrong_S,'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    subtitle('(without considering boundaries)')
    xlabel('Endowment (Y)')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','_Fig0.pdf'));

%% Graphs 2 (policy functions)

% Plot (Fig 1: Optimal consumption)
figname1 = strcat('Fig1: Optimal consumption', ' (\lambda=',num2str(lambda),')');

figure('Name',figname1)
subplot(2,2,1)
    plot(Y(tinit:tend),c1(tinit:tend),'r:',...
         Y(tinit:tend),c2(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    %titlestr = strcat('Endowment and Optimal Consumption','($\lambda$=',num2str(lambda),')');
    %title(titlestr,'interpreter','latex')
    title('Optimal Consumption')
        leg1 = strcat('$c_1$','(RRA=',num2str(gamma1),')');
        leg2 = strcat('$c_2$','(RRA=',num2str(gamma2),')');
    %legend('$c_1$ (RRA=1)', '$c_2$ (RRA=1/2)')
    legend(leg1,leg2)
    grid;

subplot(2,2,2)
    plot(Y(tinit:tend),S(tinit:tend),'k',...
         Y(tinit:tend),W1(tinit:tend),'r:',...
         Y(tinit:tend),W2(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Asset Price and Wealth')
    legend('Asset Price ($S$)','Wealth of Agent 1','Wealth of Agent 2')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','_Fig1.pdf'));    

%--------------------------------------
% Plot (Fig 2: Optimal Portfolio)
figname2 = strcat('Fig2: Optimal Portfolio', ' (\lambda=',num2str(lambda),')');

figure('Name',figname2)
subplot(2,2,1)
    plot(Y(tinit:tend), w11(tinit:tend),'r:',......
         Y(tinit:tend), w21(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Optimal Portfolio: risky asset')
    legend('$\omega_1^{(1)}$','$\omega_2^{(1)}$')
    grid;

subplot(2,2,2)
    plot(Y(tinit:tend), w12(tinit:tend),'r:',......
         Y(tinit:tend), w22(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Optimal Portfolio: riskless asset')
    legend('$\omega_1^{(2)}$','$\omega_2^{(2)}$')
    grid;

subplot(2,2,3)
    plot(Y(tinit:tend),N11(tinit:tend),'r:',...
         Y(tinit:tend),N21(tinit:tend),'b--',...
        'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    %titlestr1 = strcat('Riksy Asset Shares', '($\lambda$=',num2str(lambda),')');
    %title(titlestr1,'Interpreter','latex')
    title('Risky Asset Shares')
    legend('$N_1^{(1)}$', '$N_2^{(1)}$')
    grid;
    
subplot(2,2,4)
    plot(Y(tinit:tend),NB1(tinit:tend),'r:',...
         Y(tinit:tend),NB2(tinit:tend),'b--',...
        'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Money invested in riskless asset')
    legend('$B*N_1^{(2)}$', '$B*N_2^{(2)}$')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','_Fig2.pdf'));    

%--------------------------------------
% Plot (Fig 3: Asset Prices)
figname3 = strcat('Fig2: Asset Prices', ' (\lambda=',num2str(lambda),')');

figure('Name',figname3)
% ORIGINAL graph: 2 x 2
%{
subplot(2,2,1)
    plot(Y,S,'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    xlabel('Endowment (Y)')
    grid;

subplot(2,2,2)
    plot(Y, pd,'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Price-Dividend Ratio ($S/Y$)')
    grid;

subplot(2,2,3)
    plot(Y(tinit:tend),r(tinit:tend),'k',...
         Y(tinit:tend),-psi(tinit:tend),'r:',...
         Y(tinit:tend),beta(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Asset Prices I')
    legend('Interest Rate ($r$)', 'Price of Risk (-$\psi$)', 'Expected Rate of Return ($\beta$)')
    grid;

subplot(2,2,4)
    plot(Y(tinit:tend),m(tinit:tend),'r:',...
         Y(tinit:tend),sigmat(tinit:tend),'b--','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Asset Prices II')
    legend('Stochastic Discount Factor ($m$)', 'Stock Volatility ($\sigma_t$)')
    grid;
%}

% NEW graph: 2 x 3
subplot(2,3,1)
    plot(Y,S,'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    xlabel('Endowment (Y)')
    grid;

subplot(2,3,2)
    plot(Y, pd,'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Price-Dividend Ratio ($S/Y$)')
    grid;

subplot(2,3,3)
    plot(Y(tinit:tend),r(tinit:tend),'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Interest Rate ($r$)')
    %legend('Interest Rate ($r$)')
    grid;

subplot(2,3,4)
    plot(Y(tinit:tend),-psi(tinit:tend),'r',...
         Y(tinit:tend),beta(tinit:tend),'b','LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Asset Prices')
    legend('Price of Risk (-$\psi$)', 'Expected Rate of Return ($\beta$)')
    grid;

subplot(2,3,5)
    plot(Y(tinit:tend),m(tinit:tend),'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Stochastic Discount Factor ($m$)')
    %legend('Stochastic Discount Factor ($m$)')
    grid;

subplot(2,3,6)
    plot(Y(tinit:tend),sigmat(tinit:tend),'LineWidth', 1.5);
    xlabel('Endowment (Y)')
    title('Stock Volatility ($\sigma_t$)')
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','_Fig3.pdf'));    


%--------------------------------------
%Figure (changing "X axis")
figure('Name','Fig4')
subplot(2,2,1)
plot(s(tinit:tend), r(tinit:tend),'LineWidth', 1.5);
xlabel('Relative Consumption of Agent 2 ($s = c_2/Y$)')
title('Interest Rate ($r_t$)')
%legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;

subplot(2,2,2)
plot(s(tinit:tend), pd(tinit:tend),'LineWidth', 1.5);
xlabel('Relative Consumption of Agent 2 ($s = c_2/Y$)')
title('Price-Dividend Ratio ($S/Y$)')
%legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;

subplot(2,2,3)
plot(s(tinit:tend), sigmat(tinit:tend),'LineWidth', 1.5);
xlabel('Relative Consumption of Agent 2 ($s = c_2/Y$)')
title('Stock Volatility ($\sigma_t$)')
%legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;

subplot(2,2,4)
plot(s(tinit:tend), -psi(tinit:tend),'LineWidth', 1.5);
xlabel('Relative Consumption of Agent 2 ($s = c_2/Y$)')
title('Sharpe Ratio ($-\psi_t$)')
%legend('\omega_1^{(1)}','\omega_2^{(1)}')
grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','_Fig4.pdf'));    

%% Extra analysis
ry = (r(2:end) - r(1:end-1))/deltaY; % Derivative of "r" wrt "Y"
% by Ito's lemma
vol_r = sigma*Y(1:end-1).*ry;

figure('Name','Volatility')
plot(Y(1:end-1),vol_r, 'LineWidth',1.5)
title('Volatility of the Interest Rate')
