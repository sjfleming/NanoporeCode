function params = update_step(alpha, beta, logA, logEm, states, observations)
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
    % properly normalized so that sum(10.^gamma,1) = 1
    for t = 1:numel(observations)-1
        % xi(i,j,t) transition prob from state i to j at observation level t
        xi(:,:,t) = repmat(alpha(:,t),1,numel(states)) + logA + repmat(logEm(t+1,:)',1,numel(states)) ...
            - log10(sum(sum(10.^( repmat(alpha(:,t),1,numel(states)) + logA + repmat(logEm(t+1,:)',1,numel(states)) ))));
    end
    % properly normalized so that sum(sum(10.^(xi(:,:,t)),1)) = 1 for all t
    
    % gamma(i,t) is the probability of being in state i at observation
    % number t, given the observation sequence and the model parameters
    
    % xi(i,j,t) is the probability of making a transition from state i to
    % state j at observation number t, given the observation sequence and
    % model parameters
    
    % update the transition matrix according to baum-welch
    trans = log10(sum(10.^xi,3)) - log10(repmat(sum(sum(10.^xi,3),2),1,size(xi,2)));
    % use this to update the parameters of p_back, p_stay, p_forward, p_skip
    p_back = [];
    p_forward = [];
    p_skip = [];
    p_stay = [];
    for i = 1:size(trans,1)
        for j = 1:size(trans,2)
            switch j-i
                case 0
                    p_stay = [p_stay, 10.^trans(i,j)];
                case -1
                    p_back = [p_back, 10.^trans(i,j)];
                case 1
                    p_forward = [p_forward, 10.^trans(i,j)];
                case 2
                    p_skip = [p_skip, 10.^trans(i,j)];
            end
        end
    end
    params.p_back = mean(p_back);
    params.p_stay = mean(p_stay);
    params.p_forward = mean(p_forward);
    params.p_skip = mean(p_skip);
    
    % we want the model scaling update to reflect the full posterior, so we don't just
    % update for the most probable level, but for all level likelihoods
    % over all possible alignments. minimize the weighted residuals.
    model = repmat(arrayfun(@(x) x.level_mean, states)',1,numel(observations)); % model level means in matrix form
    data = repmat(arrayfun(@(x) x.level_mean, observations),numel(states),1); % measured level means in matrix form
    
%     lim = log10(0.2);% fit messes up if you include very unlikely points, 10% prob seems reasonable
%     x = model(gamma>lim);
%     y = data(gamma>lim);
%     e = 10.^(-gamma(gamma>lim));
%     ft = fittype('m*x + b','independent',{'x'},'coefficients',{'m','b'});
%     f = fit(x,y,ft,'StartPoint',[1, 0],'Weights',e,'Robust','bisquare');
    
    eqn = @(p) sum(sum((p(1)*model+p(2)*ones(size(model))-data).^2 .* 10.^gamma));
    options = optimset('Display', 'off');
    f = fmincon(eqn,[1;0],[-1, 0; 0, 1],[0, max(max(abs(data-model)))],[],[],[],[],[],options);
    
    params.p_scale = mean([1, f(1)]);
    params.p_offset = f(2)/2;
    
end