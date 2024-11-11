How to Build and Solve Continuous-time Heterogeneous Agents Models in Asset Pricing? The Martingale Approach and the Finite Difference Method (2024)

(1) MainCode: It solves the PDE of stock price for a two-agent model with risk-aversion heterogeneity for RRA_agent1 = 2*RRA_agent2. It implements the Finite Difference method with implicit and upwind schemes.

(2) PolicyFunctions: It uses the output of Maincode to calculate policy functions.

(3) SensitivityFunction: A function that solves the PDE of stock price and obtains policy functions. This function depends on model parameters. 

(4) SensitivityAnalysis: It uses SensitivityFunction to obtain policy functions when parameters of the endowment process and the relative risk aversion of agents change.

(5) AlternativeCalibration: It uses SensitivityFunction to solve the PDE of stock price and obtain policy functions for three different parameter calibrations.

(6) RRArange: A function that calculates the range of the RRA of each agent and the maximum value of the RRA of the more risk-averse agent when RRA_agent1 = 2*RRA_agent2.

(7) NumericalExpectation: A function that approximates the expectation presented in the wealth of the more risk-averse agent using Gaussian quadrature.

(8) PDEwithoutBoundaries: It solves the PDE of stock price without considering boundary conditions.

(9) Convergence_Stability_Monotonicity: It plot convergence, stability, and monotonicity of the numerical method.
