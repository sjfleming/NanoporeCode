function analytical_Pe(Pe,fig)

% analytical solution, Redner, Eqn 2.2.31 (also 2.3.29)
% Stephen Fleming
% 7/22/17

    % solution
    tau = @(x) (1./(2*x) - (1./(4*x.^2)) .* (1-exp(-2*x)) );
    
    % plot
    
    figure(fig)
    plot(-1*Pe,tau(Pe))

    set(gca,'yscale','log')
    set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
    ylabel('Escape time (L^2/D)')
    xlabel('-1 * Péclet number (L\sigma|V| / 2k_BT)')
    hold on

end
