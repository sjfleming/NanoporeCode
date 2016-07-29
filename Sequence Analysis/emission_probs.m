function obs_prob = emission_probs(measurement, measured_current_range, state, p_noise, p_deep)
% emission_probs(level_mean, level_stdv, stdv_mean, p_noise, p_deep)
% returns the value of the pdf for a measured current of a given level
% Inputs:
% 'measurement' the observation(s): struct that contains two fields:
%       'level_mean' level's mean current
%       'level_stdv' standard deviation of level current
% 'measured_current_range' row 2-vector of of the entire current range
% 'state' is a struct that contains three fields:
%       'level_mean' level's mean current
%       'level_stdv' standard deviation of level current
%       'stdv_mean' standard deviation of level's mean current
% 'p_noise' the a priori probability of observing a junk noise value
% 'p_deep' the a priori probability of observing a deep block level
% 'p_scale' the scale factor for the states
% 'p_offset' the offset for the states (applied after scale)
% Stephen Fleming
% 7/24/16
    
    % finds the overlap between an emission probability function and the
    % measured level.
    
    xx = linspace(measured_current_range(1), measured_current_range(2), 1000);
    pdf = zeros(size(xx));
    
    % the normal level that should be seen: use a Gaussian distribution
    p_normal = 1 - p_noise - p_deep;
    pdf = pdf + p_normal * normpdf(xx, state.level_mean, state.stdv_mean); % just the std of the mean
    
    % deep block is about 25% of the conductance of the normal level
    % use a wide Gaussian with a stdv of 10%
    pdf = pdf + p_deep * normpdf(xx, state.level_mean*0.25, state.level_mean*0.1);
    
    % noise pdf should integrate to p_noise: use a flat dist
    pdf = pdf + p_noise / diff(measured_current_range) * ...
        double(and(xx >= measured_current_range(1), xx <= measured_current_range(2))) ...
        .* ones(size(pdf));
    
    % now calculate the overlap between the emission pdf and the level
    for i = 1:numel(measurement)
        obs_prob(i) = sum( pdf .* normpdf(xx, measurement(i).level_mean, measurement(i).level_stdv) * (xx(2)-xx(1)));
    end
    
end