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
Date:   2024 (Nov)
Paper base: Wang(1996) 
----------------------------
%}
%=========================================
clear; clc;
close all;

%% Run the main m-file

run G_Main_PDE_S2F.m

iterationsV = 1:n;

%% Monotonicity

irange=1:I+1;


    % Graph
    figure('Name','Conv-Sta-Mono')
    subplot(2,3,1)
    plot(dist,'LineWidth', 1.5)
    title('Convergence')
    xlabel('Iterations (n)')
    ylabel('$||S^n - S^{n-1} ||$')
    grid;

    subplot(2,3,2)
    plot(error,'LineWidth', 1.5)
    title('Stability')
    xlabel('Iterations (n)')
    ylabel('$L_2$ Error Norm')
    grid;

    subplot(2,3,3)
    %yyaxis left
    plot(irange, -X,'b',...
         irange, -Z, 'r--','LineWidth', 1.5)
    
    %yyaxis right
    %plot(irange, -(1/deltat)*ones(I+1,1),'LineWidth', 1.5)

    title('Monotonicity')
    xlabel('Grid Points ($i$)')
    %ylabel('$-X_i$ and $-Z_i$')
    legend('$-X_i$', '$-Z_i$')
    grid;

    % Save the figure
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(gcf,'PaperPosition', [0 0 1 1]);
    print(h, '-dpdf', strcat('P1','_FigCSM.pdf'));