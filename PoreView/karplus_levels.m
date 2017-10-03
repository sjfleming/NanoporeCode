function levels = karplus_levels(data, expected_levels_per_second, false_positives_per_second, filter)
% karplus_levels(data, expected_levels_per_second, false_positives_per_second, filter)
% implements a version of the level-finding algorithm
% with details described here:
% https://gasstationwithoutpumps.wordpress.com/2014/06/17/segmenting-noisy-signals-revisited/
%
% 'data' should be two columns, the first is time, the second is current
% returns a struct 'levels' containing information about each level
% Stephen Fleming, 11/16/16

    % recursive level-finding based on a probability threshold for a
    % likelihood function
    function level_search(ind1,ind2)
        
        warning('off','MATLAB:colon:nonIntegerIndex'); % turn off warning since fminbnd uses non-integer indices
        options = optimset('TolX',1);
        [possible_transition_index, min_prob] = fminbnd(@(x) log_posterior_odds(x,ind1,ind2), ind1+minpts, ind2-minpts, options);
        
        tracker(end+1,:) = [min_prob, possible_transition_index];
        if min_prob < p
            % this index is a real transition, save it and recursively look
            % at the two new levels
            level_transition_indices(a) = possible_transition_index;
            a = a+1;
            fprintf('.');
            if possible_transition_index-ind1 > 2 * minpts
                level_search(ind1,possible_transition_index);
            end
            if ind2-possible_transition_index > 2 * minpts
                level_search(possible_transition_index,ind2);
            end
        end
        
    end
    
    % likelihood function giving probability that levels are the same,
    % conditioned on the index location of a possible transition
    function probability = log_posterior_odds(index,i1,i3)
        
        probability = (index-i1) * 0.5*log(fastvariance(i1,index)) + ...
            (i3-index) * 0.5*log(fastvariance(index,i3)) - ...
            (i3-i1) * 0.5*log(fastvariance(i1,i3)); % - ...
            %log_prior;

    end
    
    % fast variance function (I am guessing about this implementation from a hint on that blog...)
    function variance = fastvariance(ind1, ind2)
        ind2 = min(size(data,1), round(ind2));
        ind1 = max(1, round(ind1));
        num = n(ind2)-n(ind1);
        mu = (y1(ind2) - y1(ind1))/num;
        variance = (y2(ind2) - y2(ind1))/num - mu^2;
    end
    
    % precompute for fast variance calculation (make it work with nans too)
    logic = ~isnan(data(:,2));
    n = cumsum(logic);
    y1 = cumsum(data(logic,2));
    y2 = cumsum(data(logic,2).^2);
    
    % set the threshholds
    sampling = data(2,1)-data(1,1); % sampling interval (s)
    fs = 1/sampling;
    k = filter/(fs/2); % ratio of filter frequency to Nyquist frequency (if filter is Nyquist frequency, then k=1)
    %log_prior = log(expected_levels_per_second) - log(fs - expected_levels_per_second);
    %p = -log_prior + 1/k * log(false_positives_per_second/fs);
    %p = -1/k * (log(fs) - log(false_positives_per_second)); % Schreiber and Karplus 2015, 11
    p = -1/k * (log(fs - expected_levels_per_second) - log(expected_levels_per_second)); % Schreiber and Karplus 2015, 10
    %minpts = 1/k * (1/filter)/sampling; % number of data points in the shortest resolvable level with this filter setting
    minpts = (1/filter)/sampling;
    
    % run the level search (recursive)
    tracker = [];
    a = 1;
    level_transition_indices = [];
    level_search(1,size(data,1));
    
    % SJF 2017/10/03 trying to make this more reliable
    [yy,xx] = hist(tracker(:,1),-logspace(10,-3,100));
    g = fit(-log10(-xx'),yy','gauss1','startpoint',[size(tracker,1)/10,nanmean(-log10(-tracker(:,1))),nanstd(-log10(-tracker(:,1)))]);
    cutoff = -1 * 10^-norminv(0.01,g.b1,g.c1); % look at where the CDF of a normal exceeds 1%
    level_transition_indices = tracker(tracker(:,1)<cutoff,2)';
    
    % go through each level and summarize its data, keeping it in a struct
    levels = cell(numel(level_transition_indices)+1,1);
    level_transition_indices = [1, round(sort(level_transition_indices)), size(data,1)];
    for j = 1:(numel(level_transition_indices)-1)
        levels{j}.start_time = data(level_transition_indices(j),1);
        levels{j}.end_time = data(level_transition_indices(j+1),1);
        levels{j}.duration = data(level_transition_indices(j+1),1)-data(level_transition_indices(j),1);
        levels{j}.current_mean = nanmean(data(level_transition_indices(j):level_transition_indices(j+1),2));
        levels{j}.current_median = nanmedian(data(level_transition_indices(j):level_transition_indices(j+1),2));
        levels{j}.current_std = nanstd(data(level_transition_indices(j):level_transition_indices(j+1),2));
    end

    fprintf('\n')

end