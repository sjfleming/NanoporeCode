function out = viterbi_assignment(observations, states)
% viterbi_assignment(observations, states)
% assigns a sequence of observations to a given path through the possible
% state space in a hidden Markov model
% observations is a column vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a column vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% Stephen Fleming
% 7/24/16

    % state space S, numbered
    % (all possible underlying states in the hidden Markov model)
    S = (1:numel(states))'; % states, NO deep block states, NO noise state
    
    % transition matrix A, probabilities
    % (probability of going from one state to another)
    % initial guesses for transitions, out of the blue
    p_stay = max(0.1,(numel(observations)-numel(states))/numel(observations)); % a fudge for the fact that sometimes level-finding overfinds
    p_back = (1-p_stay)/12; % total probability of going back, any number of steps
    p_skip = (1-p_stay)/12; % total probability of skipping, any number of steps
    p_forward = 10*(1-p_stay)/12; % total probability of taking one forward step, what we expect to happen
    A = transition_matrix(numel(states), p_back, p_stay, p_forward, p_skip);
    
    % emission probability function, 'emission'
    % not a matrix because we don't have discrete observation states
    % (probability of an observation given the underlying state)
    p_deep = max(0.001, sum(arrayfun(@(x) x.level_means, observations)<min(states))); % a priori probability of observing a deep block
    p_noise = max(0.01, (numel(observations)-numel(states)-p_stay*numel(observations)) / numel(observations)); % a priori probability of a meaningless level in the data
    emission = @(obs,state) emission_probs(obs, [min(arrayfun(@(x) x.level_means, observations)), max(arrayfun(@(x) x.level_means, observations))], ...
        state.level_mean, state.level_stdv, state.stdv_mean, p_noise, p_deep);
    
    % initial state vector init, probabilities
    % (probabilities of starting in each state)
    init = [0.9, 0.05, 0.03, 0.012, 0.005, 0.002, 0.001, zeros(1, max(0,numel(states)-7))];
    if numel(init)>numel(states)
        init = init(1:numel(states)) / sum(init(1:numel(states))); % concatenate and fix to sum to 1
    end
    
    % initialize the big matrix T to trace path
    % rows are all possible states
    % columns are the observations
    % each element stores the log of the most likely path so far
    % and a pointer to the previous element of that path
    T = zeros(numel(states), numel(observations);
    % initial values
    T(:,1) = init.*emission(observations(1));
    
end