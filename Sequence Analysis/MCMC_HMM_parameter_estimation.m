function out = MCMC_HMM_parameter_estimation(observations, states, param_start, proposalmethod, T, N)
% MCMC_HMM_parameter_estimation(observations, states, param_start, proposalmethod, T, N)
% samples from the posteriors for the parameters of an HMM using a
% step-size T (equivalent to temperature in a simulated-annealing sense).
% observations is a vector struct of observed current levels
% observations struct contains 'level_mean' and 'level_stdv' fields
% states is a vector struct of model current levels
% states struct contains 'level_mean', 'level_stdv', and 'stdv_mean' fields
% param_start struct contains 'p_stay', 'p_back', 'p_skip', 'p_forward',
% 'p_deep', and 'p_noise', or can be an empty set [].
% proposalmethod is a string, either 'classic', 'baumwelch', 'traceback',
% or 'classicbaumwelch'
% T is a step-size (simulated annealing temp).  Note that we sample from
% the true posterior since our step size depends on T but the acceptance
% probability is for T=1.
% N is the desired number of samples.

% Idea for suboptimal path sampling (here used to generate the MCMC "step")
% was taken from Sean R. Eddy, "Multiple alignment using hidden Markov
% models," ISMB 1995.

% Stephen Fleming
% 8/1/16
    
    % initializing parameters of the HMM
    if isempty(param_start)
        p.p_stay = mean(cellfun(@(x) max(0.1,(numel(x)-numel(states))/numel(x)), observations)); % a fudge for the fact that sometimes level-finding overfinds
        p.p_back = (1-p.p_stay)/12; % total probability of going back, any number of steps
        p.p_skip = (1-p.p_stay)/12; % total probability of skipping, any number of steps
        p.p_forward = 10*(1-p.p_stay)/12; % total probability of taking one forward step, what we expect to happen
        p.p_deep = mean(cellfun(@(y) max(0.001, sum(arrayfun(@(x) x.level_mean, y) < 0.8*min(arrayfun(@(x) x.level_mean, states))) / numel(y)), observations)); % a priori probability of observing a deep block
        p.p_noise = mean(cellfun(@(x) max(0.01, (numel(x)-numel(states)-p.p_stay*numel(x)) / numel(x)), observations)); % a priori probability of a meaningless level in the data
        p.scale = 1;
        p.offset = 0;
    else
        try
            p.p_stay = param_start.p_stay;
            p.p_back = param_start.p_back;
            p.p_skip = param_start.p_skip;
            p.p_forward = param_start.p_forward;
            p.p_noise = param_start.p_noise;
            p.p_deep = param_start.p_deep;
            p.scale = param_start.scale;
            p.offset = param_start.offset;
        catch ex
            display('Parameters incorrectly assigned.')
        end
    end
    
    samples{1} = p;
    
    % initial sample
    % calculate initial state probs, transition matrix, and emission matrix
    A = transition_matrix(numel(states), p.p_back, p.p_stay, p.p_forward, p.p_skip);
    logA = log10(A);
    init = [0.9, 0.05, 0.03, 0.012, 0.005, 0.002, 0.001, 1e-10*ones(1, max(0,numel(states)-7))];
    if numel(init)>numel(states)
        init = init(1:numel(states)) / sum(init(1:numel(states))); % concatenate and fix to sum to 1
    end
    logInit = log10(init);

    % generate the exponentiated forward variable
    for k = 1:numel(observations)
        I_range = [min(arrayfun(@(x) x.level_mean, observations{k})), max(arrayfun(@(x) x.level_mean, observations{k}))];
        emission = @(obs,state) emission_probs(obs, I_range, state, p.p_noise, p.p_deep);
        lem = log10(cell2mat(arrayfun(@(y) arrayfun(@(x) emission(x,y), observations{k}), states, 'uniformoutput', false)')');
        logEm{k} = lem;
        
        [nlf, nlb, ~, nla, nlc] = forward_backward_variables_scaled(observations{k}, states, logInit, logA, logEm{k}, T);
        logF{k} = nlf;
        logbeta{k} = nlb;
        logalpha{k} = nla;
        logC{k} = nlc;
    end
    
    % keep track of acceptance ratio
    a = 0;
    r = 0;
    
    % proposals for variables
    function new_params = classicProposal(p, stdv, iter, var)
        % generate a classical proposal by using a normal distribution
        % around one of the current parameters and renormalizing others
        np(1) = p.p_forward;
        np(2) = p.p_back;
        np(3) = p.p_stay;
        np(4) = p.p_skip;
        
        % one variable random step
%         whichvar = mod(iter,4)+1; % integer 1 through 4
        whichvar = var;
        np(whichvar) = np(whichvar) + randn(1)*stdv;
        
        % all variables random step
%         np = np + randn(1,4)*stdv;
        
        % two variable tradeoff
%         whichvars = randsample(1:4,2); % two non-repeated integers 1 through 4
%         step = randn(1)*stdv;
%         np(whichvars(1)) = np(whichvars(1)) + step;
%         np(whichvars(2)) = np(whichvars(2)) - step;
        
        % random change in random number of variables
%         whichvars = randsample(1:4,randi(4)); % random number of non-repeated integers 1 through 4
%         np(whichvars) = np(whichvars) + randn(1,numel(whichvars))*stdv;
        
        % make sure it's within bounds
        if any(np<=0)
            np = 1-np;
            np = np / sum(np);
            np = 1-np;
        end
        np = np / sum(np);
        
        new_params.p_forward = np(1);
        new_params.p_back = np(2);
        new_params.p_stay = np(3);
        new_params.p_skip = np(4);
    end
    
    function new_params = baumwelchProposal(p, iter, var)
        % generate a classical proposal by using a normal distribution
        % around one of the current parameters and renormalizing others
        np(1) = p.p_forward;
        np(2) = p.p_back;
        np(3) = p.p_stay;
        np(4) = p.p_skip;
%         whichvar = mod(iter,4)+1; % integer 1 through 4
        
        temp_logA = log10(transition_matrix(numel(states), p.p_back, p.p_stay, p.p_forward, p.p_skip));
        % baum welch update for that variable
        ra = randi(numel(observations),1); % use a random observation to generate this sample
        [alpha_bar, beta_bar, ~, ~, ~] = forward_backward_variables_scaled(observations{ra}, states, logInit, temp_logA, logEm{ra}, 1);
        bwparam = update_step(alpha_bar, beta_bar, temp_logA, logEm{ra}, states, observations{ra}, 0);
        switch var
            case 1
                np(1) = bwparam.p_forward;
            case 2
                np(2) = bwparam.p_back;
            case 3
                np(3) = bwparam.p_stay;
            case 4
                np(4) = bwparam.p_skip;
        end
        
        % make sure it's within bounds
        np = np / sum(np);
        new_params.p_forward = np(1);
        new_params.p_back = np(2);
        new_params.p_stay = np(3);
        new_params.p_skip = np(4);
    end
    
    function new_params = tracebackProposal(p, T, iter)
        % generate a backtrace proposal by using a backtrace through one of
        % the molecules at a given temperature and calculating one of the
        % parameters based on this
        np(1) = p.p_forward;
        np(2) = p.p_back;
        np(3) = p.p_stay;
        np(4) = p.p_skip;
        %whichvar = mod(iter,4)+1; % integer 1 through 4
        
        % probabilistic traceback: get a new sample path
        % initialization
        ra = randi(numel(observations),1); % use a random observation to generate this sample
        state_inds = zeros(1,numel(observations{ra}));
        state_probs = logF{ra}(:,end); % sums to one, already normalized by the scaling procedure in 'forward_backward_variables_scaled'
        state_inds(end) = numel(states) - (find(rand(1)>=[flipud(cumsum(10.^state_probs)); 0],1,'first') - 1) + 1;
        % induction
        for t = numel(observations{ra}):-1:2
            state_probs = logF{ra}(:,t-1) + (1/T) * logA(:,state_inds(t)) - log10(sum(10.^(logF{ra}(:,t-1) + (1/T) * logA(:,state_inds(t)))));
            state_inds(t-1) = numel(states) - (find(rand(1)>=[flipud(cumsum(10.^state_probs)); 0],1,'first') - 1) + 1;
        end
        
        %switch whichvar
            %case 1
                np(1) = max(1e-5,sum(diff(state_inds)==1)) /(numel(state_inds)-1);
            %case 2
                np(2) = max(1e-5,sum(diff(state_inds)<0))  /(numel(state_inds)-1);
            %case 3
                np(3) = max(1e-5,sum(diff(state_inds)==0)) /(numel(state_inds)-1);
            %case 4
                np(4) = max(1e-5,sum(diff(state_inds)>1))  /(numel(state_inds)-1);
        %end
        
        np = np / sum(np);
        new_params.p_forward = np(1);
        new_params.p_back = np(2);
        new_params.p_stay = np(3);
        new_params.p_skip = np(4);
    end
    
    % iterative sample generation via MCMC
    for iter = 2:N+1
        
        % get a proposal sample
        switch proposalmethod
            case 'classic'
                new_params = classicProposal(p, 0.06, iter);
            case 'traceback'
                new_params = tracebackProposal(p, T, iter);
            case 'baumwelch'
                new_params = baumwelchProposal(p, iter);
            case 'classicbaumwelch'
                va = mod(iter,4)+1;
                new_params = classicProposal(p, 0.3, iter, va);
                new_params = baumwelchProposal(new_params, iter, va);
        end
        
        % calculate probability of observations given the new parameters
        % just for a select few
        new_logA = log10(transition_matrix(numel(states), new_params.p_back, new_params.p_stay, new_params.p_forward, new_params.p_skip));
        testers = randsample(1:numel(observations),3); % pick 3 to look at
        % for the current sample
        for k = testers
            [nlf, nlb, ~, nla, nlc] = forward_backward_variables_scaled(observations{k}, states, logInit, logA, logEm{k}, T);
            logF{k} = nlf;
            logbeta{k} = nlb;
            logalpha{k} = nla;
            logC{k} = nlc;
        end
        % for the proposed sample
        new_logF = logF;
        new_logbeta = logbeta;
        new_logalpha = logalpha;
        new_logC = logC;
        for k = testers
            [nlf, nlb, ~, nla, nlc] = forward_backward_variables_scaled(observations{k}, states, logInit, new_logA, logEm{k}, T);
            new_logF{k} = nlf;
            new_logbeta{k} = nlb;
            new_logalpha{k} = nla;
            new_logC{k} = nlc;
        end
        
        % calculate the acceptance probability
        Prob_Ratio = 10.^( sum(cellfun(@(x) sum(x), new_logC(testers))) - sum(cellfun(@(x) sum(x), logC(testers))) );
        if isnan(Prob_Ratio)
            pause
        end
        
        % accept or reject move
        if rand(1)<=Prob_Ratio
            % accept!
            
            % update the current point to the proposed sample
            logF = new_logF;
            logC = new_logC;
            logA = new_logA;
            logbeta = new_logbeta;
            logalpha = new_logalpha;
            p = new_params;
            
            % draw stuff
%             figure(1)
%             clf
%             image(logF{ra},'cdatamapping','scaled')
%             hold on
%             plot(1:numel(observations{ra}),state_inds,'r-','LineWidth',3)
%             drawnow;
            figure(2)
            clf
            plot(cellfun(@(x) x.p_forward, samples),'g')
            hold on
            plot(cellfun(@(x) x.p_back, samples),'r')
            plot(cellfun(@(x) x.p_stay, samples),'k')
            plot(cellfun(@(x) x.p_skip, samples),'c')
            ylim([0 1])
            drawnow;
            a = a+1;
            display(['Accepted: prob_ratio ' num2str(Prob_Ratio) ', AR = ' num2str(round(a/(a+r)*100)) '%'])
            
        else
            r = r+1;
            display(['Rejected: prob_ratio ' num2str(Prob_Ratio) ', AR = ' num2str(round(a/(a+r)*100)) '%'])
        end
        
        % take a sample, always!
        samples{iter} = p;
        
    end
    
    out = samples;
    
end