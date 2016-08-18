function [observations, hidden_states, hidden_state_ind, noise_ind, deep_ind] = generate_HMM_observations(states, p)
% generate_HMM_observations(states, p) creates a series of HMM
% states and their corresponding observations based on input probabilities
% of the HMM.  used for HMM parameter estimation testing purposes.
% states is a struct witht the following fields:
%   pA is a vector of mean current levels
%   pAstd is a vector of the current std of those levels
%   std_mean is a vector of the std of many measurements of a mean level
% p is a struct containing the HMM probabilities and parameters:
%   p_back
%   p_stay
%   p_forward
%   p_skip
%   p_noise
%   p_deep
%   scale
%   offset

% Stephen Fleming
% 8/1/16
    
    % inputs
    pA = arrayfun(@(x) x.level_mean, states);
    pAstd = arrayfun(@(x) x.level_stdv, states);
    std_mean = arrayfun(@(x) x.stdv_mean, states);
    
    % hidden model states
    i = 1;
    state_inds = 1;
    while i<=numel(pA)
        n = rand(1);
        if n < p.p_back && i>=2
            state_inds(end+1) = i-1;
            i = i-1;
        elseif n < p.p_back+p.p_stay
            state_inds(end+1) = i;
        elseif n < p.p_back+p.p_stay+p.p_forward
            state_inds(end+1) = i+1;
            i = i+1;
        else % skip
            state_inds(end+1) = i+2;
            i = i+2;
        end
    end
    state_inds = state_inds(state_inds<=numel(pA));
    hidden_state_ind = state_inds;
    hidden_states = struct('level_mean',num2cell(pA(state_inds)),'level_stdv',num2cell(pAstd(state_inds)),'stdv_mean',num2cell(std_mean(state_inds)));
    
    % observations
    noise_ind = rand([1,numel(state_inds)]) < p.p_noise;
    deep_ind = rand([1,numel(state_inds)]) < p.p_deep;
    noise_ind = noise_ind & ~deep_ind;
    levs = pA(state_inds) + randn(size(pA(state_inds))) .* std_mean(state_inds);
    levs(noise_ind) = rand([1,sum(noise_ind)])*(max(pA)-min(pA)) + min(pA);
    levs(deep_ind) = levs(deep_ind)*0.25;
    observations = struct('level_mean',num2cell(levs * p.scale + p.offset),'level_stdv',num2cell(pAstd(state_inds) * p.scale));
    
end