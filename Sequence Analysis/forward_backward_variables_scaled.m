function [alpha_bar_T, beta_bar_T, c_T, alpha_bar, c] = forward_backward_variables_scaled(observations, states, logInit, logA, logEm, T)
% forward_backward_variables_scaled(observations, states, logInit, logA, logEm)
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
% returns the probability of observing 'observation' from the model state
% 'state' (NOT log prob)

% Stephen Fleming
% 8/1/16
    
    % alpha columns are observations
    % alpha rows are states
    alpha_bar_T = zeros(numel(states),numel(observations));
    alpha_bar = zeros(numel(states),numel(observations));
    
    % initialization
    alpha_bar_T(:,1) = (1/T) * (logInit + logEm(1,:));
    alpha_bar(:,1) = logInit + logEm(1,:);
    c_T = zeros(1,numel(observations)); % c is the scale factors, also in log space
    c = zeros(1,numel(observations));
    
    % induction
    for t = 1:(numel(observations)-1) % step through observations
        
        for j = 1:numel(states) % for each state
            
            alpha_bar_T(j,t+1) = log10(sum( 10.^(alpha_bar_T(:,t) + (1/T) * logA(:,j)) )) + (1/T) * logEm(t+1,j); % log space
            alpha_bar(j,t+1) = log10(sum( 10.^(alpha_bar(:,t) + logA(:,j)) )) + logEm(t+1,j); % log space
            
        end
        
        % rescale, and keep track of scaling
        c_T(t+1) = log10(sum(10.^alpha_bar_T(:,t+1)));
        alpha_bar_T(:,t+1) = alpha_bar_T(:,t+1) - c_T(t+1);
        c(t+1) = log10(sum(10.^alpha_bar(:,t+1)));
        alpha_bar(:,t+1) = alpha_bar(:,t+1) - c(t+1);
        
    end
    
    % beta columns are observations
    % beta rows are states
    beta_bar_T = zeros(numel(states),numel(observations));
    
    % initialization
    beta_bar_T(:,numel(observations)) = zeros(numel(states),1); % since log10(1)=0
    
    % induction
    for t = (numel(observations)-1):-1:1 % step through observations, backward
        
        for i = 1:numel(states) % for each state
            
            beta_bar_T(i,t) = log10(sum( 10.^(beta_bar_T(:,t+1) + (1/T) * logA(i,:)' + (1/T) * logEm(t+1,:)') )) - c_T(t); % log space, rescaled according to alpha_bar scaling
            
        end
        
    end
    
end