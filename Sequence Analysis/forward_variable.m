function alpha = forward_variable(observations, states, logInit, logA, logEm)
% forward_variable(observations, states, init, A, emission)
% calculates the value of the forward variable
% for a given hidden Markov model
% see Rabiner 1989, IEEE 77(2), Section III.A

% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% logInit is the log10 initial probabilities
% logA is the log10 transition matrix probabilities
% logEm is the log10 precomputed emission matrix
% returns the probability of observing 'observation' from the model state
% 'state' (NOT log prob)

% Stephen Fleming
% 7/28/16
    
    % alpha columns are observations
    % alpha rows are states
    alpha = zeros(numel(states),numel(observations));
    
    % initialization
    alpha(:,1) = logInit + logEm(1,:);
    
    % induction
    for t = 1:(numel(observations)-1) % step through observations
        
        for j = 1:numel(states) % for each state
            
            alpha(j,t+1) = log10(sum( 10.^(alpha(:,t) + logA(:,j)) )) + logEm(t+1,j); % log space
            
        end
        
    end
    
    % output
    % alpha gets returned
    
end