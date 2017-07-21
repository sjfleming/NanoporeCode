% Stephen Fleming
% 7/21/17
% Gardiner Eqn. 5.2.160
% mean first passage time
% one reflecting boundary (0), one absorbing (L)

%% stuff

% inputs
V = 0.03; % voltage in Volts
s = 0.5; % how many electronic charges per base
T = 273.15 + 25; % K
n = 14; % number of bases past constriction

% constants and derived quantities
e_0 = 1.6e-19; % coulombs
nu = 0.59; % Flory exponent
kB = 1.38e-23; % SI units
Lp = 1e-9; % Tinland 1997
Lb = 0.5e-9; % length per base ssDNA, meters
L = n * Lb; % length of domain, meters
sigma = s * e_0 / Lb; % charge per unit length

%% quasipotential

% electrostatic
U_e = @(x) sigma * V * x;

% entropic
U_s = @(x) nu * kB * T * (log((x + Lp)/L) + log(1 - x/(L + Lp)));

% total
U = @(x) U_e(x) + U_s(x);

figure(1)
xx = linspace(0,L,1000);
plot(xx * 1e9,U_e(xx) / (kB*T),'k')
hold on
plot(xx * 1e9,U_s(xx) / (kB*T),'k--')
plot(xx * 1e9,U(xx) / (kB*T),'k','linewidth',3)
xlim([0, Inf])
xlabel('Domain (nm)')
ylabel('Quasipotential (kT)')

%% force

% electrostatic
F_e = -1 * sigma * V;

% entropic
F_s = @(x) -1 * nu * kB * T * (1./((x + Lp)/L)*(1/L) + 1./(1 - x/(L + Lp))*(-1/(L + Lp)));

% total
F = @(x) F_e + F_s(x);

figure(2)
plot(xx * 1e9,F_e*ones(size(xx)) / (kB*T / Lb),'k')
hold on
plot(xx * 1e9,F_s(xx) / (kB*T / Lb),'k--')
plot(xx * 1e9,F(xx) / (kB*T / Lb),'k','linewidth',3)
xlim([0, Inf])
xlabel('Domain (nm)')
ylabel('Electrokinetic force (kT / L_b)')

%% probability distribution for initial position

w = @(x) (sigma*0.16/(kB*T)) * exp(-1*sigma*0.16/(kB*T) * x); % normalized, 160mV equilibrium dist.

figure(3)
plot(xx * 1e9,w(xx),'k')
xlabel('Domain (nm)')
ylabel('Probability density')
title('Starting location')

%% integrating factor

psi = @(x) exp(2/(kB*T) * integral(@(y) F(y),0,x)); % use matlab to integrate numerically
