function out = viterbi_assignment(observations, states)
% viterbi_assignment(observations, states)
% assigns a sequence of observations to a given path through the possible
% state space in a hidden Markov model
% observations is a column vector of observed current levels
% states is a column vector of model current levels

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
    figure()
    image(A,'cdatamapping','scaled')
    title('Transition probability matrix')
    
    % emission probability function E
    % not a matrix because we don't have discrete observation states
    % (probability of an observation given the underlying state)
    p_deep = max(0.001, sum(observations<min(states))); % a priori probability of observing a deep block
    p_noise = max(0.01, (numel(observations)-numel(states)-p_stay*numel(observations)) / numel(observations)); % a priori probability of a meaningless level in the data
    
    
    % initial state vector init, probabilities
    % (probabilities of starting in each state)
    
    
    % initialize the big matrices T1 and T2 to trace path
    % rows are all possible states
    % columns are the observations
    
    
end