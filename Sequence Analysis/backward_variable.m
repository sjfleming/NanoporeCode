function beta = backward_variable(observations, states, logA, logEm)
% backward_variable(observations, states, init, A, emission)
% calculates the value of the forward variable
% for a given hidden Markov model
% see Rabiner 1989, IEEE 77(2), Section III.A

% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% logA is the log10 transition matrix probabilities
% logEm is the log10 precomputed emission matrix
% returns the probability of observing 'observation' from the model state
% 'state' (NOT log prob)

% Stephen Fleming
% 7/27/16
    
    % beta columns are observations
    % beta rows are states
    beta = zeros(numel(states),numel(observations));
    
    % initialization
    beta(:,numel(observations)) = zeros(numel(states),1); % since log10(1)=0
    
    % induction
    for t = (numel(observations)-1):-1:1 % step through observations, backward
        
        for i = 1:numel(states) % for each state
            
            beta(i,t) = log10(sum( 10.^(beta(:,t+1) + logA(i,:)' + logEm(t+1,:)') )); % log space
            
        end
        
    end
    
    % output
    % beta gets returned
    
end