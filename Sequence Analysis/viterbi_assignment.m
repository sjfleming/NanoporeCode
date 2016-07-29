function out = viterbi_assignment(observations, states)
% viterbi_assignment(observations, states)
% assigns a sequence of observations to a given path through the possible
% state space in a hidden Markov model
% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields

% Stephen Fleming
% 7/24/16
    
    % transition matrix A, probabilities
    % (probability of going from one state to another)
    % initial guesses for transitions, out of the blue
    p_stay = max(0.1,(numel(observations)-numel(states))/numel(observations)); % a fudge for the fact that sometimes level-finding overfinds
    p_back = (1-p_stay)/12; % total probability of going back, any number of steps
    p_skip = (1-p_stay)/12; % total probability of skipping, any number of steps
    p_forward = 10*(1-p_stay)/12; % total probability of taking one forward step, what we expect to happen
    A = transition_matrix(numel(states), p_back, p_stay, p_forward, p_skip);
    logA = log10(A);
    
    % emission probability function, 'emission'
    % not a matrix because we don't have discrete observation states
    % (probability of an observation given the underlying state)
    p_deep = max(0.001, sum(arrayfun(@(x) x.level_mean, observations) < 0.8*min(arrayfun(@(x) x.level_mean, states))) / numel(observations)); % a priori probability of observing a deep block
    p_noise = max(0.01, (numel(observations)-numel(states)-p_stay*numel(observations)) / numel(observations)); % a priori probability of a meaningless level in the data
    I_range = [min(arrayfun(@(x) x.level_mean, observations)), max(arrayfun(@(x) x.level_mean, observations))];
    p_scale = 1;
    p_offset = 0;
    emission = @(obs,state) emission_probs(obs, I_range, state, p_noise, p_deep);
    % pre-compute all values for speed
    logEm = log10(cell2mat(arrayfun(@(y) arrayfun(@(x) emission(x,y), observations), states, 'uniformoutput', false)')');
    
    % initial state vector init, probabilities
    % (probabilities of starting in each state)
    init = [0.9, 0.05, 0.03, 0.012, 0.005, 0.002, 0.001, 1e-10*ones(1, max(0,numel(states)-7))];
    if numel(init)>numel(states)
        init = init(1:numel(states)) / sum(init(1:numel(states))); % concatenate and fix to sum to 1
    end
    
    % initialize the big matrix T to trace path
    % rows are all possible states
    % columns are the observations
    % dimension 3 contains: prob, pointer i, pointer j in that order
    T = zeros(numel(states), numel(observations), 3);
    % initial values
    T(:,1,1) = log10( init .* arrayfun(@(x) emission(observations(1), x), states) );
    
    % fill the big matrix
    % each element (i,j,1) stores the log prob of the most likely path so far
    % and a pointer (i,j,2:3) to the previous element of that path
    for j = 2:numel(observations) % columns are observations
        
        for i = 1:numel(states) % rows are all possible states
            
            % all possible previous states times transition to this state
            % times emission for this state
            [m, ind] = max( T(:,j-1,1) + logA(:,i) + logEm(j,i) ); % sum log probabilities
            T(i,j,1) = m; % the probability of the maximally probable path to here
            T(i,j,2) = ind; % row index
            T(i,j,3) = j-1; % column index, always the previous column...
            
        end
        
    end
    
    % get the most probable path by finding the last most probable state
    [~,ind] = max(T(:,end,1)); % state index
    z = nan(1,numel(observations)); % best path state indices
    z(end) = ind;
    
    % trace back through the big matrix to get the sequence of states
    for j = numel(observations):-1:2
        
        z(j-1) = T(z(j),j,2); % pointer to the previous row index, i.e. state index
        
    end
    
    % figure out which states are deep and which are noise
    % since the emission probability subsumes deep blocks and noise in this
    % model, we need to go back and see if each level was normal, deep
    % block, or noise
    state_type = cell(1,numel(observations));
    for i = 1:numel(z)
        reg = emission_probs(observations(i), I_range, states(z(i)), 0, 0);
        noi = emission_probs(observations(i), I_range, states(z(i)), 0.5, 0);
        dee = emission_probs(observations(i), I_range, states(z(i)), 0, 0.5);
        if and(reg>noi, reg>dee)
            state_type{i} = 'normal';
        elseif and(noi>reg, noi>dee)
            state_type{i} = 'noise';
        elseif and(dee>reg, dee>noi)
            state_type{i} = 'deep';
        else
            state_type{i} = 'undetermined';
        end
    end
    
    % package the output
    out.state_indices = z;
    out.state_sequence = arrayfun(@(x) states(x), z);
    out.state_type = state_type;
    out.log_prob_matrix = T;
    out.observations = observations;
    
end