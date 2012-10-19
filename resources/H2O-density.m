function H2Odensity
% Fitting procedure based on
% http://www.eng.cam.ac.uk/help/tpl/programs/Matlab/curve_fitting.html
% Kevin L. Davies, 4/9/11 (modified 11/22/11)

m = 0.01801528; % Specific mass / (kg/mol) [McBride2002, NIST2009]
T_degC = 0:85; % Temperature in deg C
T = T_degC + 273.15; % Temperature / K
n_Takenaka = @(T)(1 - (T - 277.13152).^2.*(T + 123.03534).*(T - 240.86147)./(609628.6*(T - 190.02667).*(T - 242.90545)))*999.9734/m; % Density / (mol/m3)
% From Takenaka, M. & Masui, R. "Measurement of the Thermal Expansion of Pure Water in the Temperature Range 0 degC to 85 degC", Metrologia, 1990, 27, 165

n = n_Takenaka(T);
bestcoeffs = fminsearch(@fun, [1 1 1 1], [], T, n);
a = bestcoeffs(1)
b = bestcoeffs(2)
c = bestcoeffs(3)
d = bestcoeffs(4)
n_fun = @(T)a*(b*T + c).^3 + d; % Density / (mol/m3)
n_fit = n_fun(T);
%n_fun2 = @(T_degC)999.9734 - 31.3625/209.199*(T_degC - 4).^2; % Density / (mol/m3)
%n_fit2 = n_fun2(T_degC);

% Now compare n with n_fit.
plot(T, n, '-', T, n_fit, '.');
%hold on
%plot(T, n_fit2, 'r.');
legend('Takenaka & Masui 1990','Fit')
%legend('Takenaka & Masui 1990', 'Fit1', 'Fit2')
title('Density of Liquid H_2O')
xlabel('Temperature / K')
ylabel('Density / mol m^-3')

function out = fun(coeff, T, n)
a = coeff(1);
b = coeff(2);
c = coeff(3);
d = coeff(4);
n_fun = @(T)a*(b*T + c).^3 + d; % Density / (mol/m3)
DIFF = n_fun(T) - n;
SQ_DIFF = DIFF.^2;

out = sum(SQ_DIFF);
