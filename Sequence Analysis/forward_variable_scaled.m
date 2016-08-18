function [alpha_bar, c] = forward_variable_scaled(observations, states, logInit, logA, logEm, T)
% forward_variable_scaled(observations, states, logInit, logA, logEm, T)
% calculates the value of the forward variable
% for a given hidden Markov model
% see Rabiner 1989, IEEE 77(2), Sections III.A and V.A

% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% logInit is the log10 initial probabilities
% logA is the log10 transition matrix probabilities
% logEm is the log10 precomputed emission matrix
% returns the log10 probability of observing 'observation' from the model
% state 'state'
% T is a temperature by which alpha_bar is exponentiated alpha^(1/T)
% T=1 gives the normal forward variable

% Stephen Fleming
% 8/1/16
    
    % alpha columns are observations
    % alpha rows are states
    alpha_bar = zeros(numel(states),numel(observations));
    
    % initialization
    alpha_bar(:,1) = (1/T) * (logInit + logEm(1,:));
    c = zeros(1,numel(observations)); % c is the scale factors, also in log space
    
    % induction
    for t = 1:(numel(observations)-1) % step through observations
        
        for j = 1:numel(states) % for each state
            
            alpha_bar(j,t+1) = log10(sum( 10.^(alpha_bar(:,t) + (1/T) * logA(:,j)) )) + (1/T) * logEm(t+1,j); % log space
            
        end
        
        % rescale every time, and keep track of scaling
        s = log10(sum(10.^alpha_bar(:,t+1)));
        c(t+1) = s;
        alpha_bar(:,t+1) = alpha_bar(:,t+1) - s;
        
    end
    
    % output
    % alpha gets returned
    
end