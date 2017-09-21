function [kappa, pulses] = pulse_correlation(sigdata,tr,chan,levels)
% pulse_correlation
% Stephen Fleming
% 9/21/17
    
    % pulse timings
    pulses = [];
    d = sigdata.getViewData(tr);
    logic = diff(medfilt1(d(:,4),2))>0.5;
    logic = logic & ~[false; logic(1:end-1)]; % one point every spike
    inds = find(logic); % indices of those spikes
    inds = inds([true; diff(inds)>10]);
    dt = d(2,1)-d(1,1);
    for i = 1:numel(inds)
        ind = (tr(1) + inds(i)*dt)/sigdata.si; % index of raw data
        pulses(end+1) = sigdata.si * sigdata.findNext(@(d) d(:,chan)>0.05,ind-100);
    end
    
    % cumulative distribution function for level durations
    kappa = correlation_metric(pulses, cellfun(@(x) x.start_time, levels), mean(cellfun(@(x) x.duration, levels)));
    
end