function kappa = correlation_metric(pulse_timings, level_change_timings, tau)
% a metric for correlation between pulses and level changes
% Stephen Fleming
% 2/26/16
    
    % cumulative distribution function for level durations
    function c = cdf(t,tau)
        c = 1 - exp(-1*t/tau);
    end
    
    % calculate p value for each pulse
    p = nan(1,numel(pulse_timings));
    for i = 1:numel(pulse_timings)-1
        ind = find(level_change_timings>pulse_timings(i)&level_change_timings<pulse_timings(i+1),1,'first'); % closest level change before next pulse
        if ~isempty(ind)
            p(i) = cdf(level_change_timings(ind)-pulse_timings(i), tau);
        end
    end
    
    % calculate kappa
    bins = ceil(numel(pulse_timings)/10);
    bincenters = 1/bins/2 : 1/bins : 1 - 1/bins/2;
    n = hist(p, bincenters);
    kappa = sqrt( std(n)^2*bins / numel(p)^2 );
    
    % plot if no arguments passed out
    if nargout==0
        figure(1)
        clf(1)
        bar(bincenters,n)
        xlabel('p value')
        ylabel('Number of pulses')
        title(['kappa = ' num2str(kappa,2)])
        set(gca,'fontsize',20)
        display(kappa)
    end
    
end