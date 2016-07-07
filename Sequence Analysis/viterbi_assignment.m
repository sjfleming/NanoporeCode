function out = viterbi_assignment(observations, states)
% viterbi_assignment(observations, states)
% assigns a sequence of observations to a given path through the possible
% state space in a hidden Markov model
% observations is a column vector of observed current levels
% states is a column vector of model current levels

    % state space S, numbered
    % (all possible underlying states in the hidden Markov model)
    S = (1:(2*numel(states)+1))'; % states, deep block states, a noise state
    
    % transition matrix A, probabilities
    % (probability of going from one state to another)
    p_noise = 0.01; % a priori probability of an observation with no meaning
    p_deep = 0.02; % a priori probability of transitioning to a deep block
    p_back = 0.1; % total probability of going back, any number of steps
    p_skip = 0.1; % total probability of skipping, any number of steps
    p_stay = 0.05; % a fudge for the fact that sometimes level-finding overfinds
    p_forward = 1 - p_noise - p_deep - p_back - p_skip - p_stay; % probability of doing what should actually happen
    A = zeros(numel(S),numel(S));
    block = zeros(numel(states),numel(states));
    for i = 1:numel(states)
        row = zeros(1,numel(states));
        % stay
        row(i) = p_stay;
        if i<numel(states)
            row(i+1) = p_forward;
        end
        % skip
        skipinds = [i+2, numel(states)];
        if skipinds(1)<=numel(states)
            skipinds = min(numel(states),skipinds);
            skipinds = skipinds(1):skipinds(2);
            skip = (p_skip.^(1:numel(skipinds)))/sum(p_skip.^(1:numel(skipinds)))*p_skip;
            row(skipinds) = skip;
        end
        % back
        backinds = [1, i-1];
        if backinds(2)>=1
            backinds = max(1,backinds);
            backinds = backinds(1):backinds(2);
            back = (fliplr(p_back.^(1:numel(backinds))))/sum(p_back.^(1:numel(backinds)))*p_back;
            row(backinds) = back;
        end
        % noise
        
        % deep
        
        % ensure A is normalized even though it should be close
        block(i,1:numel(states)) = row/(sum(row));
    end
    % A(i,j) is the probability of a transition from state i to state j
    A(1:numel(states),1:numel(states)) = block * (1-p_noise-p_deep);
    A(numel(states)+1:2*numel(states),1:numel(states)) = block * (1-p_noise-p_deep); % deeps to regulars
    % deep
    A(1:numel(states),numel(states)+1:2*numel(states)) = block * p_deep;
    A(numel(states)+1:2*numel(states),numel(states)+1:2*numel(states)) = block * p_deep;
    % noise
    A(:,numel(S)) = p_noise;
    % leave the noise state ambiguous... will have to deal with noise->other later
    
    % emission matrix E, probabilities
    % (probability of an observation given the underlying state)
    a;
    
    % initial state vector init, probabilities
    % (probabilities of starting in each state)
    
    
    % initialize the big matrices T1 and T2 to trace path
    % rows are all possible states
    % columns are the observations
    
    
end