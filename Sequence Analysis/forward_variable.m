function alpha = forward_variable(observations, states, init, A, emission)
% forward_variable(observations, states, init, A, emission)
% calculates the value of the forward variable
% for a given hidden Markov model
% see Rabiner 1989, IEEE 77(2), Section III.A

% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% init is the initial probabilities
% A is the transition matrix
% emission is a function handle and emission(observation, state)
% returns the probability of observing 'observation' from the model state
% 'state'

% Stephen Fleming
% 7/26/16
    
    % alpha columns are observations
    % alpha rows are states
    
    % initialization
    alpha(:,1) = init .* arrayfun(@(x) emission(observations, x), states);
    
    % induction
    
    
    % output
    
    
end