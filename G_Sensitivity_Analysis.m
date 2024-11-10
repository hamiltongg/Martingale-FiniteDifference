%% P1
% Author: Hamilton Galindo
% Started (Mar 2021) - Updated (Feb,Mar,April 2023, July 2024)
% Two-agent model
% Martingale approach
%==========================================================================

clear;
close all;

p = mfilename('fullpath');           % the full path of this m-file
[filepath,name,ext] = fileparts(p);  % Split the "full path": We need "name" (of this m-file) 

%% Parameters
% [C1]mu - (1/2)sigma^2 > 0 : it makes E(dy/y)>0 (+ growth) : see Eq.2, Wang(1996)
% [C2]growth condition:
    % rho > (1/2)*max(0,mu - (1/4)sigma^2)
    rho_vec = [0.02 0.1];

for i=1:2    
rho     = rho_vec(i);  %impatience rate ==> discount factor = e^(-rho*t) 
%mu      = 0.05; %E[dY/Y]: baseline model
%sigma   = 0.3;  %Volatility term of dY/Y: baseline model

lambda = 2/3;
b = 4*( (1-lambda)^2 )/lambda^2; % coming from the agent's 1 budget constraint at t=0

%% Limits of mu (given rho & using C1 and C2)
sigmavec = 0.1:0.1:1;
lower_mu = (1/2)*sigmavec.^2; 
%upper_mu = rho + (1/4)*sigmavec.^2;
upper_mu = 2*rho + (1/4)*sigmavec.^2;

subplot(2,2,i)
plot(sigmavec,lower_mu,...
    sigmavec,upper_mu,...
    'LineWidth',1.5)
legend('$lower_{\mu} = (1/2)\sigma^2$','$upper_{\mu} = 2\rho + (1/4)\sigma^2$','interpreter','latex','Location','best')
    titlestr = strcat('Limits of $\mu$ ','($\rho$=',num2str(rho),')');
    title(titlestr,'interpreter','latex')
xlabel('$\sigma$')
ylabel('$\mu$')
grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','Sen_Fig0.pdf'));    
end

%% Analysis 1: Fix \mu = 0.05
rho    = rho_vec(2); % I choose \rho=0.1 bc it gives more room for \mu and \sigma
lambda = 2/3;

mu = 0.05;
sigma_vec = [0.1 0.2];
gamma1 = 1;

for i=1:2
    result_mufix{i} = G_function_Main_PDE_S2F(rho, lambda, mu,sigma_vec(i),gamma1);
end

%% Analysis 2: Fix \sigma = 0.2
rho    = rho_vec(2); % I choose \rho=0.1 bc it gives more room for \mu and \sigma
lambda = 2/3;

mu_vec = [0.05 0.1];
sigma = 0.2;
gamma1 = 1;

for i=1:2
    result_sigmafix{i} = G_function_Main_PDE_S2F(rho, lambda, mu_vec(i),sigma,gamma1);
end

%% Analysis 3: \gamma1 
rho    = rho_vec(2); % I choose \rho=0.1 bc it gives more room for \mu and \sigma
lambda = 2/3;

mu = 0.05;
sigma = 0.2;
gamma1_vec = [1, 2.2];

for i=1:2
    result_gamma1{i} = G_function_Main_PDE_S2F(rho, lambda, mu,sigma,gamma1_vec(i));
end

%% Graphs: Fix \mu = 0.05

Y = result_mufix{1}.Y;
tinit = 2;
tend = 500;

% Plot (Fig 3: Asset Prices)
figname3 = strcat('Fig2: Asset Prices', ' ($\mu=$',num2str(result_mufix{1}.mu),')');

figure('Name',figname3)
subplot(2,3,1)
    plot(Y(tinit:tend),result_mufix{1}.vars{1,3}(tinit:tend,1), ...
         Y(tinit:tend),result_mufix{2}.vars{1,3}(tinit:tend,1), 'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\sigma =$',num2str(result_mufix{1}.sigma));
        leg2 = strcat('$\sigma =$',num2str(result_mufix{2}.sigma));
    legend(leg1,leg2)
    grid;

subplot(2,3,2)
    plot(Y(tinit:tend),result_mufix{1}.vars{1,3}(tinit:tend,2), ...
         Y(tinit:tend),result_mufix{2}.vars{1,3}(tinit:tend,2), 'LineWidth',1.5)
    title('Price-Dividend Ratio ($S/Y$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\sigma =$',num2str(result_mufix{1}.sigma));
        leg2 = strcat('$\sigma =$',num2str(result_mufix{2}.sigma));
    legend(leg1,leg2)
    grid;

subplot(2,3,3)
    plot(Y(tinit:tend),result_mufix{1}.vars{1,3}(tinit:tend,3), ...
         Y(tinit:tend),result_mufix{2}.vars{1,3}(tinit:tend,3), 'LineWidth',1.5)
    title('Interest Rate ($r$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\sigma =$',num2str(result_mufix{1}.sigma));
        leg2 = strcat('$\sigma =$',num2str(result_mufix{2}.sigma));
    legend(leg1,leg2)
    grid;

subplot(2,3,4)
    plot(Y(tinit:tend),-result_mufix{1}.vars{1,3}(tinit:tend,4), ...
         Y(tinit:tend),-result_mufix{2}.vars{1,3}(tinit:tend,4), 'LineWidth',1.5)
    title('Price of Risk (-$\psi$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\sigma =$',num2str(result_mufix{1}.sigma));
        leg2 = strcat('$\sigma =$',num2str(result_mufix{2}.sigma));
    legend(leg1,leg2)
    grid;

subplot(2,3,5)
    plot(Y(tinit:tend),result_mufix{1}.vars{1,3}(tinit:tend,5), ...
         Y(tinit:tend),result_mufix{2}.vars{1,3}(tinit:tend,5), 'LineWidth',1.5)
    title('Expected Rate of Return ($\beta$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\sigma =$',num2str(result_mufix{1}.sigma));
        leg2 = strcat('$\sigma =$',num2str(result_mufix{2}.sigma));
    legend(leg1,leg2)
    grid;

subplot(2,3,6)
    %Y,result_mufix{1}.vars{1,3}(:,6),'r', ...
    %Y,result_mufix{2}.vars{1,3}(:,6),'r:',...
    plot(Y(tinit:tend),result_mufix{1}.vars{1,3}(tinit:tend,7), ...
         Y(tinit:tend),result_mufix{2}.vars{1,3}(tinit:tend,7),...
         'LineWidth',1.5)
    title('Stock Volatility ($\sigma_t$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\sigma =$',num2str(result_mufix{1}.sigma));
        leg2 = strcat('$\sigma =$',num2str(result_mufix{2}.sigma));
    legend(leg1,leg2)
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','Sen_Fig1.pdf'));    

%% Graphs: Fix \sigma = 0.2

% Plot (Fig 3: Asset Prices)
figname4 = strcat('Fig3: Asset Prices', ' ($\sigma=$',num2str(result_sigmafix{1}.sigma),')');

figure('Name',figname4)
subplot(2,3,1)
    plot(Y(tinit:tend),result_sigmafix{1}.vars{1,3}(tinit:tend,1), ...
         Y(tinit:tend),result_sigmafix{2}.vars{1,3}(tinit:tend,1), 'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\mu =$',num2str(result_sigmafix{1}.mu));
        leg2 = strcat('$\mu =$',num2str(result_sigmafix{2}.mu));
    legend(leg1,leg2)
    grid;

subplot(2,3,2)
    plot(Y(tinit:tend),result_sigmafix{1}.vars{1,3}(tinit:tend,2), ...
         Y(tinit:tend),result_sigmafix{2}.vars{1,3}(tinit:tend,2), 'LineWidth',1.5)
    title('Price-Dividend Ratio ($S/Y$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\mu =$',num2str(result_sigmafix{1}.mu));
        leg2 = strcat('$\mu =$',num2str(result_sigmafix{2}.mu));
    legend(leg1,leg2)
    grid;

subplot(2,3,3)
    plot(Y(tinit:tend),result_sigmafix{1}.vars{1,3}(tinit:tend,3), ...
         Y(tinit:tend),result_sigmafix{2}.vars{1,3}(tinit:tend,3), 'LineWidth',1.5)
    title('Interest Rate ($r$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\mu =$',num2str(result_sigmafix{1}.mu));
        leg2 = strcat('$\mu =$',num2str(result_sigmafix{2}.mu));
    legend(leg1,leg2)
    grid;

subplot(2,3,4)
    plot(Y(tinit:tend),-result_sigmafix{1}.vars{1,3}(tinit:tend,4), ...
         Y(tinit:tend),-result_mufix{2}.vars{1,3}(tinit:tend,4), 'LineWidth',1.5)
    title('Price of Risk (-$\psi$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\mu =$',num2str(result_sigmafix{1}.mu));
        leg2 = strcat('$\mu =$',num2str(result_sigmafix{2}.mu));
    legend(leg1,leg2)
    grid;

subplot(2,3,5)
    plot(Y(tinit:tend),result_sigmafix{1}.vars{1,3}(tinit:tend,5), ...
         Y(tinit:tend),result_sigmafix{2}.vars{1,3}(tinit:tend,5), 'LineWidth',1.5)
    title('Expected Rate of Return ($\beta$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\mu =$',num2str(result_sigmafix{1}.mu));
        leg2 = strcat('$\mu =$',num2str(result_sigmafix{2}.mu));
    legend(leg1,leg2)
    grid;

subplot(2,3,6)
    %Y,result_sigmafix{1}.vars{1,3}(:,6),'r', ...
    %Y,result_sigmafix{2}.vars{1,3}(:,6),'r:',...
    plot(Y(tinit:tend),result_sigmafix{1}.vars{1,3}(tinit:tend,7), ...
         Y(tinit:tend),result_sigmafix{2}.vars{1,3}(tinit:tend,7),...
         'LineWidth',1.5)
    title('Stock Volatility ($\sigma_t$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\mu =$',num2str(result_sigmafix{1}.mu));
        leg2 = strcat('$\mu =$',num2str(result_sigmafix{2}.mu));
    legend(leg1,leg2)
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','Sen_Fig2.pdf'));    

%% Graphs: \gamma1

% Plot (Fig 3: Asset Prices)
figname4 = strcat('Fig3: Asset Prices',...
           ' ($\sigma=$',num2str(result_sigmafix{1}.sigma),')',...
           ' ($\mu=$',num2str(result_mufix{1}.mu),')');

figure('Name',figname4)
subplot(2,3,1)
    plot(Y(tinit:tend),result_gamma1{1}.vars{1,3}(tinit:tend,1), ...
         Y(tinit:tend),result_gamma1{2}.vars{1,3}(tinit:tend,1), 'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\gamma_1 =$',num2str(result_gamma1{1}.gamma1));
        leg2 = strcat('$\gamma_1 =$',num2str(result_gamma1{2}.gamma1));
    legend(leg1,leg2)
    grid;

subplot(2,3,2)
    plot(Y(tinit:tend),result_gamma1{1}.vars{1,3}(tinit:tend,2), ...
         Y(tinit:tend),result_gamma1{2}.vars{1,3}(tinit:tend,2), 'LineWidth',1.5)
    title('Price-Dividend Ratio ($S/Y$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\gamma_1 =$',num2str(result_gamma1{1}.gamma1));
        leg2 = strcat('$\gamma_1 =$',num2str(result_gamma1{2}.gamma1));
    legend(leg1,leg2)
    grid;

subplot(2,3,3)
    plot(Y(tinit:tend),result_gamma1{1}.vars{1,3}(tinit:tend,3), ...
         Y(tinit:tend),result_gamma1{2}.vars{1,3}(tinit:tend,3), 'LineWidth',1.5)
    title('Interest Rate ($r$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\gamma_1 =$',num2str(result_gamma1{1}.gamma1));
        leg2 = strcat('$\gamma_1 =$',num2str(result_gamma1{2}.gamma1));
    legend(leg1,leg2)
    grid;

subplot(2,3,4)
    plot(Y(tinit:tend),-result_gamma1{1}.vars{1,3}(tinit:tend,4), ...
         Y(tinit:tend),-result_gamma1{2}.vars{1,3}(tinit:tend,4), 'LineWidth',1.5)
    title('Price of Risk (-$\psi$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\gamma_1 =$',num2str(result_gamma1{1}.gamma1));
        leg2 = strcat('$\gamma_1 =$',num2str(result_gamma1{2}.gamma1));
    legend(leg1,leg2)
    grid;

subplot(2,3,5)
    plot(Y(tinit:tend),result_gamma1{1}.vars{1,3}(tinit:tend,5), ...
         Y(tinit:tend),result_gamma1{2}.vars{1,3}(tinit:tend,5), 'LineWidth',1.5)
    title('Expected Rate of Return ($\beta$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\gamma_1 =$',num2str(result_gamma1{1}.gamma1));
        leg2 = strcat('$\gamma_1 =$',num2str(result_gamma1{2}.gamma1));
    legend(leg1,leg2)
    grid;

subplot(2,3,6)
    %Y,result_gamma1{1}.vars{1,3}(:,6),'r', ...
    %Y,result_gamma1{2}.vars{1,3}(:,6),'r:',...
    plot(Y(tinit:tend),result_gamma1{1}.vars{1,3}(tinit:tend,7), ...
         Y(tinit:tend),result_gamma1{2}.vars{1,3}(tinit:tend,7),...
         'LineWidth',1.5)
    title('Stock Volatility ($\sigma_t$)')
    xlabel('Endowment (Y)')
        leg1 = strcat('$\gamma_1 =$',num2str(result_gamma1{1}.gamma1));
        leg2 = strcat('$\gamma_1 =$',num2str(result_gamma1{2}.gamma1));
    legend(leg1,leg2)
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('Wang1996','Sen_Fig3.pdf'));    