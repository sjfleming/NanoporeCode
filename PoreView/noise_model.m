function y = noise_model(x,temp,conductance,current,filter,Ra,Cm,loss_tangent,Cin,V_headstage)
% noise_model is a full workup of nanopore noise sources
% Stephen Fleming
% 7/28/17

% typical input values:
% temp = 22; % in Celsius
% conductance = 2; % in nS
% Ra = 3e7; % access resistance = rho/4*a (J.E. Hall, 1975), rho = 0.0895ohm*m for 1M KCl
% Cm = 0.45e-12; % typical membrane capacitance, in Farads
% loss_tangent = 1; % Dave ref to http://pubs.acs.org/doi/pdf/10.1021/jp0680138
% current = 1e-9; % in amperes
% Cin = 4e-12; % headstage input capacitance is about 4pF (Sakmann and Neher say 15pF, p.112)
% V_headstage = 3e-9; % input voltage noise on headstage op-amp = 3nV/sqrt(Hz)

    R = 1/(conductance*1e-9);
    k = 1.38 * 10^-23; % Boltzmann constant
    T = temp + 273.15; % absolute temperature in Kelvin
    Rf = 500e6; % feedback resistor in Axopatch with beta = 1 in whole cell mode
    dielectric_loss = @(f) 8*pi*k*T*f*Cm*loss_tangent * 1e18;
    ion = 1.12e-9 * (current*1e-12)^2 * (1e9)^2; % ion number fluctuations proportional to current squared (10/20/16 data for MspA from Oxford)
    % Y = @(f,Ra,Cm) sqrt(-1)*2*pi*f*Cin + 1./( 2*Ra + 1./( 1/R + sqrt(-1)*2*pi*f.*Cm ) ) + 1./Rf; % complex conductance (admittance)
    Y = @(f,Ra,Cm) 1./( 2*Ra + 1./( 1/R + sqrt(-1)*2*pi*f.*Cm ) ) + 1./Rf; % complex conductance (admittance)
    % headstage = V_headstage^2.*real(Y(f,Ra,Cm).*conj(Y(f,Ra,Cm))) * 10^18;
    headstage = @(f) V_headstage^2.*(2*pi*f*Cin).^2 * 10^18; % Sakmann-Neher p.112
    % the filter
    cutoff = filter;
    factor = pi/1.8031;
    [b,a] = besself(4,cutoff*factor);
    h = freqs(b,a,x);
    noise_model = 4*k*T*real(Y(x,Ra,Cm)*1e18.*h.*conj(h)) + ion*real(h.*conj(h)) + headstage(x).*real(h.*conj(h)) + dielectric_loss(x).*real(h.*conj(h));
    
    y = noise_model;

end