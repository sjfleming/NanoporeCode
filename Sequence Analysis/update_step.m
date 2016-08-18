function params = update_step(alpha, beta, logA, logEm, states, observations, scaleoffsetflag)
% update_step(alpha, beta)
% returns updated model parameters based on the forward and backward
% probabilites of the HMM.
% 'alpha' is the log10 forward variable matrix
% 'beta' is the log10 backward variable matrix

% Stephen Fleming
% 7/28/16

    % the transition matrix only changes according to updated parameters
    % p_back, p_stay, p_forward, and p_skip
    
    % the emission matrix only changes according to updated parameters
    % p_noise, p_deep, p_scale, and p_offset
    
    % so we don't update these things directly
    % we calculate updates to the parameters
    % and use the parameters to recalculate transition and emission probs
    
    % our notion of emission probability is different, continuous
    % we want the emission probability in each state to be maximal
    % basically just fit p_scale and p_offset to the residuals
    % using the probability of a observations being in model states
    
    gamma = (alpha + beta) - repmat(log10( sum(10.^(alpha + beta), 1) ), size(alpha,1), 1);
    
    % gamma(i,t) is the probability of being in state i at observation
    % number t, given the observation sequence and the model parameters
    
    % update the transition matrix according to baum-welch
    for i = 1:numel(states)
        normalization = log10(sum(sum(10.^(repmat(alpha(i,1:end-1),numel(states),1)+repmat(logA(i,:),numel(observations)-1,1)'+logEm(2:end,:)'+beta(:,2:end)))));
        for j = 1:numel(states) % (Rabinov, eqn 95)
            trans(i,j) = log10(sum(10.^(alpha(i,1:end-1)+logA(i,j)+logEm(2:end,j)'+beta(j,2:end)))) - normalization;
        end
    end
    % use this to update the parameters of p_back, p_stay, p_forward, p_skip
    p_back = [];
    p_forward = [];
    p_skip = [];
    p_stay = [];
    for i = 4:size(trans,1)-3
        p_stay(end+1) = 10.^trans(i,i);
        p_forward(end+1) = 10.^trans(i,i+1);
        p_skip(end+1) = sum(10.^trans(i,i+2:end));
        p_back(end+1) = sum(10.^trans(i,1:i-1));
    end
    c = sum([mean(p_back), mean(p_stay), mean(p_forward), mean(p_skip)]);
    params.p_back = mean(p_back)/c;
    params.p_stay = mean(p_stay)/c;
    params.p_forward = mean(p_forward)/c;
    params.p_skip = mean(p_skip)/c;
    
    if scaleoffsetflag

        % we want the model scaling update to reflect the full posterior, so we don't just
        % update for the most probable level, but for all level likelihoods
        % over all possible alignments. minimize the weighted residuals.
        model = repmat(arrayfun(@(x) x.level_mean, states)',1,numel(observations)); % model level means in matrix form
        data = repmat(arrayfun(@(x) x.level_mean, observations),numel(states),1); % measured level means in matrix form

        eqn = @(p) sum(sum(abs(p(1)*model+p(2)*ones(size(model))-data) .* 10.^gamma));
        options = optimset('Display', 'off');
        f = fmincon(eqn,[1;0],[-1, 0; 0, 1],[0, max(max(abs(data-model)))],[],[],[],[],[],options);

        params.scale = f(1);%mean([1, f(1)]);
        params.offset = f(2);%/2;

    end
    
end