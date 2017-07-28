function [xx, tau] = first_passage_time(V,s,d,t,n,numpts)

    % inputs
    %V = 0.05; % voltage in Volts
    %s = 0.5; % how many electronic charges per base
    %d = 0.2; % diffusion constant in units of D_0
    T = 273.15 + t; % K
    %n = 14; % number of bases past constriction

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
    
    % forces
    F_e = -1 * sigma * V;
    F_s = @(x) -1 * nu * kB * T * (1./((x + Lp)/L)*(1/L) + 1./(1 - x/(L + Lp))*(-1/(L + Lp)));
    F = @(x) F_e + F_s(x);
    
    % evaluate
    xx = linspace(0,L,numpts);
    dx = xx(2) - xx(1);
    F2 = F(xx);
    psi = arrayfun(@(z) exp(2/(kB*T) * sum(F2(1:find(xx<=z,1,'last'))) * dx ), xx);
    cumsumpsi = cumsum(psi) * dx;
    tau = arrayfun(@(k) 2/D * sum(cumsumpsi(find(xx<=k,1,'last'):end) ./ psi(find(xx<=k,1,'last'):end) ) * dx, xx);

end