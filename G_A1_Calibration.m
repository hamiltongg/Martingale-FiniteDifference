%% P1
% Author: Hamilton Galindo
% Started (Mar 2021) - Updated (Feb,Mar,April 2023, July 2024)
% Two-agent model
% Martingale approach
%==========================================================================
clear;
close all;

%% Calibrations

names = {'Longstaff and Wang (2012)','Schneider (2022)', 'Chan and Kogan (2002)'};
         % col1: longstaff2012 | col2: schneider2022 | col3: chan2002
lambda = [1/2 , 2/3    , 2/3   ];
rho    = [0.01, 0.001/4, 0.0521];
mu     = [0.03, 0.0055 , 0.018 ];
sigma  = [0.1 , 0.019  , 0.0402];
gamma1 = [2   , 10     , 3     ];

%% Conditions

% Condition: mu > # (Eq. 76)
muHigherThan = 0.5*gamma1.*sigma.^2 - rho./(gamma1 - 1);
% Condition: gamma1 < Max
for i=1:3
    condition{i} = gammaRange(rho(i),mu(i),sigma(i));
end

%% Model Solution
for i=1:3
    result_cal{i} = G_function_Main_PDE_S2F(rho(i), lambda(i), mu(i), sigma(i), gamma1(i));
end

%% Graph
Y = result_mufix{1}.Y;
tinit = 2;
tend = 500;

figure('Name','Asset Prices')
subplot(2,2,1)
    plot(Y(tinit:tend),result_cal{1}.vars{1,3}(tinit:tend,1), ...
         Y(tinit:tend),result_cal{2}.vars{1,3}(tinit:tend,1),...
         Y(tinit:tend),result_cal{3}.vars{1,3}(tinit:tend,1),'LineWidth',1.5)
    title('Risky Asset Price ($S$)')
    xlabel('Endowment (Y)')
    legend(names{1},names{2},names{3})
    grid;

subplot(2,2,2)
    plot(Y(tinit:tend),result_cal{1}.vars{1,3}(tinit:tend,2), ...
         Y(tinit:tend),result_cal{2}.vars{1,3}(tinit:tend,2),...
         Y(tinit:tend),result_cal{3}.vars{1,3}(tinit:tend,2),'LineWidth',1.5)
    title('Price-Dividend Ratio ($S/Y$)')
    xlabel('Endowment (Y)')
    legend(names{1},names{2},names{3})
    grid;

% Save the figure
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
print(h, '-dpdf', strcat('P1','G_A1_Calibration.pdf'));  