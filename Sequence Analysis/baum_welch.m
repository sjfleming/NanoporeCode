function out = baum_welch(observations, states, params)
% baum_welch(observations, states, params)
% maximizes the likelihood of any possible observation sequence in a
% state space in a hidden Markov model
% by iteratively updating the parameters of the HMM.
% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields
% params struct contains 'p_stay', 'p_back', 'p_skip', 'p_forward',
% 'p_deep', and 'p_noise', or can be an empty set [].

% Stephen Fleming
% 7/28/16
    
    % parameters of the HMM
    if isempty(params)
        p_stay = max(0.1,(numel(observations)-numel(states))/numel(observations)); % a fudge for the fact that sometimes level-finding overfinds
        p_back = (1-p_stay)/12; % total probability of going back, any number of steps
        p_skip = (1-p_stay)/12; % total probability of skipping, any number of steps
        p_forward = 10*(1-p_stay)/12; % total probability of taking one forward step, what we expect to happen
        p_deep = max(0.001, sum(arrayfun(@(x) x.level_mean, observations) < 0.8*min(arrayfun(@(x) x.level_mean, states))) / numel(observations)); % a priori probability of observing a deep block
        p_noise = max(0.01, (numel(observations)-numel(states)-p_stay*numel(observations)) / numel(observations)); % a priori probability of a meaningless level in the data
        scale = 1;
        offset = 0;
    else
        try
            p_stay = params.p_stay;
            p_back = params.p_back;
            p_skip = params.p_skip;
            p_forward = params.p_forward;
            p_noise = params.p_noise;
            p_deep = params.p_deep;
            scale = params.scale;
            offset = params.offset;
        catch ex
            display('Parameters incorrectly assigned.')
        end
    end
    
    % transition matrix A, probabilities
    % (probability of going from one state to another)
    % initial guesses for transitions, out of the blue
    A = transition_matrix(numel(states), p_back, p_stay, p_forward, p_skip);
    logA = log10(A);
    
    % emission probability function, 'emission'
    % not a matrix because we don't have discrete observation states
    % (probability of an observation given the underlying state)
    I_range = [min(arrayfun(@(x) x.level_mean, observations)), max(arrayfun(@(x) x.level_mean, observations))];
    emission = @(obs,state) emission_probs(obs, I_range, state, p_noise, p_deep);
    % pre-compute all values for speed
    logEm = log10(cell2mat(arrayfun(@(y) arrayfun(@(x) emission(x,y), observations), states, 'uniformoutput', false)')');
    
    % initial state vector init, probabilities
    % (probabilities of starting in each state)
    init = [0.9, 0.05, 0.03, 0.012, 0.005, 0.002, 0.001, 1e-10*ones(1, max(0,numel(states)-7))];
    if numel(init)>numel(states)
        init = init(1:numel(states)) / sum(init(1:numel(states))); % concatenate and fix to sum to 1
    end
    logInit = log10(init);
    
    iterate = true;
    tolerance = 1e-3;
    
    ps = p_stay;
    pf = p_forward;
    pb = p_back;
    psk = p_skip;
    scaleoffsetflag = true;
    offset = 10;
    
    figure(1)
    clf
    hold on
    plot(arrayfun(@(x) x.level_mean, observations),'r','LineWidth',5)
    
    while iterate
        
        % calculate the scaled forward and backward variables
        [alpha_bar, beta_bar, c] = forward_backward_variables_scaled(observations, states, logInit, logA, logEm, 1);
        
        % update parameters!
        params = update_step(alpha_bar, beta_bar, logA, logEm, states, observations, scaleoffsetflag);
        
        % for now don't bother changing p_noise or p_deep
        % should i let these vary for individual levels???  might lose control
        
        % updates
        if scaleoffsetflag
            for i = 1:numel(states)
                states(i).level_mean = states(i).level_mean * params.scale + params.offset;
                states(i).level_stdv = states(i).level_stdv * params.scale;
                states(i).stdv_mean = states(i).stdv_mean * params.scale;
            end
        end
        
        logEm = log10(cell2mat(arrayfun(@(y) arrayfun(@(x) emission(x,y), observations), states, 'uniformoutput', false)')');
        logA = log10(transition_matrix(numel(states), params.p_back, params.p_stay, params.p_forward, params.p_skip));
        
        ps = [ps, params.p_stay];
        pf = [pf, params.p_forward];
        pb = [pb, params.p_back];
        psk = [psk, params.p_skip];
        if scaleoffsetflag
            offset = [offset, params.offset];
            if abs(offset(end)-offset(end-1)) < 0.1 && numel(offset)>5
                scaleoffsetflag = false;
            end
        end
        display(params)
        
        vit = viterbi_assignment(observations, states);
        figure(1)
        plot(arrayfun(@(x) x.level_mean, vit.state_sequence))
        drawnow;
        figure(2)
%         image(10.^logA,'cdatamapping','scaled')
        plot(ps,'c')
        hold on
        plot(pf,'m')
        plot(pb,'b')
        plot(psk,'k')
        ylim([0 1])
        drawnow;
        
        if numel(ps)>20%p(end)-p(end-1) < tolerance || 
            iterate = false;
        end
        
    end
    
    out = vit;
    out.params = params;
    
end