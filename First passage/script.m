% Stephen Fleming
% 7/21/17
% Gardiner Eqn. 5.2.160
% mean first passage time
% one reflecting boundary (0), one absorbing (L)

%% stuff

% inputs
V = 0.05; % voltage in Volts
s = 0.62; % how many electronic charges per base
d = 0.16; % diffusion constant in units of D_0
T = 273.15 + 25; % K
n = 14; % number of bases past constriction
tolerance = 1e-3; % for numerical integration, accurate to this

% constants and derived quantities
e_0 = 1.6e-19; % coulombs
nu = 0.59; % Flory exponent
D_0 = 3.31e-10; % diffusion const of 1.5nm sphere, m^2/s
D = d * D_0;
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
ylabel('Quasipotential (k_BT)')
set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

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
set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

%% probability distribution for initial position

w = @(x) (sigma*0.16/(kB*T)) * exp(-1*sigma*0.16/(kB*T) * x); % normalized, 160mV equilibrium dist.

figure(3)
plot(xx * 1e9,w(xx),'k')
xlabel('Domain (nm)')
ylabel('Probability density')
xlim([0, Inf])
title('Starting location')
set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

%% integrating factor

psi = @(x) exp(2/(kB*T) * arrayfun(@(z) integral(@(y) F(y),0,z,'RelTol',tolerance,'AbsTol',tolerance), x)); % use matlab to integrate numerically

%% first passage time

helperintegral = @(y) arrayfun(@(k) integral(@(z) psi(z), 0, k,'RelTol',tolerance,'AbsTol',tolerance), y);

tau = @(a) 2/D * arrayfun(@(k) integral(@(y) 1./psi(y) .* helperintegral(y), k, L,'RelTol',tolerance,'AbsTol',tolerance), a); % numerical integration

figure(4)
x1 = [0 3 6 8 10 11 12 12.5 13 13.5 13.75 14]*Lb;
plot(x1*1e9,tau(x1)*1000,'ok-')
xlabel('Starting position (nm)')
ylabel('Escape time (ms)')
set(gca,'fontsize',14,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])

%% mean first passage time

mean_tau = integral(@(a) tau(a).*w(a), 0, L,'RelTol',tolerance,'AbsTol',tolerance);

%% we need to do this more coarsely, it takes forever with this numerical integration

% just sum things
numpts = 100000;
xx = linspace(0,L,numpts);
dx = xx(2) - xx(1);
F2 = F(xx);
psi2 = arrayfun(@(z) exp(2/(kB*T) * sum(F2(1:find(xx<=z,1,'last'))) * dx ), xx);
helperintegral2 = cumsum(psi2) * dx;
tau2 = arrayfun(@(k) 2/D * sum(helperintegral2(find(xx<=k,1,'last'):end) ./ psi2(find(xx<=k,1,'last'):end) ) * dx, xx);
