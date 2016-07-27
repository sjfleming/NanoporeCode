function alpha = forward_variable(observations, states, init, A, emission)
% forward_variable(observations, states, init, A, emission)
% calculates the value of the forward variable
% for a given hidden Markov model
% see Rabiner 1989, IEEE 77(2), Section III.A

% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% init is the initial probabilities (NOT log prob)
% A is the transition matrix probabilities (NOT log prob)
% emission is a function handle and emission(observation, state)
% returns the probability of observing 'observation' from the model state
% 'state' (NOT log prob)

% Stephen Fleming
% 7/26/16
    
    % alpha columns are observations
    % alpha rows are states
    alpha = zeros(numel(states),numel(observations));
    
    % initialization
    alpha(:,1) = log10( init .* arrayfun(@(x) emission(observations(1), x), states) );
    
    % induction
    for t = 1:(numel(observations)-1) % step through observations
        
        for j = 1:numel(states) % for each state
            
            alpha(j,t+1) = log10(sum( 10.^(alpha(:,t) + log10(A(:,j))) )) + log10(emission(observations(t+1), states(j))); % log space
            
        end
        
    end
    
    % output
    % alpha gets returned
    
end