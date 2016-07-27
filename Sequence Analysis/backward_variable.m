function beta = backward_variable(observations, states, A, emission)
% backward_variable(observations, states, init, A, emission)
% calculates the value of the forward variable
% for a given hidden Markov model
% see Rabiner 1989, IEEE 77(2), Section III.A

% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% A is the transition matrix probabilities (NOT log prob)
% emission is a function handle and emission(observation, state)
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
        
        em = log10(arrayfun(@(x) emission(observations(t+1),x), states)'); % emission prob of next observation for each state
        
        for i = 1:numel(states) % for each state
            
            beta(i,t) = log10(sum( 10.^(beta(:,t+1) + log10(A(i,:)') + em) )); % log space
            
        end
        
    end
    
    % output
    % beta gets returned
    
end