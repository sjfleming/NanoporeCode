function analytical(V,s,d,t,n,fig)

% analytical solution, Redner, Eqn 2.2.31 (also 2.3.29)
% Stephen Fleming
% 7/22/17

    % inputs
    %s = 0.62; % how many electronic charges per base
    %d = 0.16; % diffusion constant in units of D_0

    T = 273.15 + t; % K
    % n = 14; % number of bases past constriction

    % constants and derived quantities
    e_0 = 1.6e-19; % coulombs
    nu = 0.59; % Flory exponent
    D_0 = 3.31e-10; % diffusion const of 1.5nm sphere, m^2/s
    D = d * D_0;
    kB = 1.38e-23; % SI units
    Lb = 0.5e-9; % length per base ssDNA, meters
    L = n * Lb; % length of domain, meters
    sigma = s * e_0 / Lb; % charge per unit length

    % solution
    
    tau = @(x) L^2/D * (-kB*T./(sigma*L*x) - (kB*T./(sigma*L*x)).^2 .* (1-exp(sigma*L*x/(kB*T))) );
    %V = linspace(0.030,0.070,500);
    
    % plot
    
    figure(fig)
    plot(V*1000,tau(V)*1000)

    set(gca,'yscale','log')
    set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
    ylabel('Escape time (ms)')
    xlabel('Voltage (mV)')
    ylim([5e-2 5e3])
    xlim([25 75])
    hold on

end
