function discreteData = find_discrete_levels(sigdata, trange, finalFrequency, medianRepeats, cutoff, pValue)
% discreteData = find_discrete_levels(trange, cutoff) uses
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.  Starts from an assumption of one level and adds
%   levels based on the size of the derivative of severely median-filtered
%   data.
%   cutoff = 0.1 works pretty well usually
%   Stephen Fleming, August 3, 2015
    
    %% divide data up into sections of 100 seconds maximum, filter heavily and downsample
    time = [];
    current = [];
    for i = 1:ceil(diff(trange)/100)
        
        tr = [trange(1)+(i-1)*100, min(trange(2),trange(1)+i*100)];
        display(['Working on time range [' num2str(tr(1)) ', ' num2str(tr(2)) ']'])
        data = sigdata.getByTime(tr);
        t = data(:,1);
        c = data(:,2)*1000;
        d = [t, c];
        clear data c t
        
        lp = filt_lpb(d,4,1000);
        [filtered] = filt_severe_median_downsample(lp,medianRepeats,(1/sigdata.si)/finalFrequency);
        t = filtered(:,1);
        c = filtered(:,2);
        clear d filtered
        
        time = [time; t];
        current = [current; c];
        
    end
    
    % get a low-pass only version of the data, not median filtered
    lp_chan = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,finalFrequency),['Low pass Bessel filter, ' num2str(finalFrequency) 'Hz']);
    lp = sigdata.getViewData(trange);
    if sum(isnan(lp(:,lp_chan(1))))>10 % if there is no view data because there are so few points
        lp = sigdata.getByTime(trange); % load the real points
    end
    lp = [lp(:,1), lp(:,lp_chan(1))];
    
%     % Add an offset manually
%     m = 0.0065; % slope of offset, pA / second
%     t_init = trange(1); % seconds, where offset correction starts
%     offset_lp = m*(lp(:,1)-t_init) * 1e-3;
%     lp(:,2) = lp(:,2) + offset_lp;
%     offset = m*(time-t_init);
%     current = current + offset;
    
    %% find level changes all together

    % take derivative
    smoothing = 10;
    deriv = smooth_derivative(current, smoothing);

    % start up the arrays of currents and times
    level_currents = mean(current);
    level_timing = [time(1) time(end)];
    discrete_current = ones(numel(time),1)*level_currents(1);

    % find peaks in the derivative and sort them
    [deriv_peaks, deriv_peak_locs] = findpeaks(abs(deriv),'MinPeakDist',2*smoothing,'MinPeakHeight',0.25*std(deriv));
    [~, sort_ind] = sort(deriv_peaks,'descend');
    num = 1;
    metric = 50;

    while and(metric(end)>cutoff, num<=numel(sort_ind))
    %for i = 1:numel(deriv_peaks)

        % find where level would split
        split_time = time(deriv_peak_locs(sort_ind(num)));
        split_level_ind = find(split_time>level_timing(:,1),1,'last');

        % check to make sure levels are statistically different
%         ts_ind1 = find(time>level_timing(split_level_ind,1),1,'first');
%         tf_ind1 = find(time<level_timing(split_level_ind,2),1,'last');
%         [~,p,~] = kstest2(current(ts_ind1:deriv_peak_locs(sort_ind(num))),current(deriv_peak_locs(sort_ind(num)):tf_ind1));
        ts_ind1 = find(lp(:,1)>level_timing(split_level_ind,1),1,'first');
        tm_ind1 = find(lp(:,1)>split_time,1,'first');
        tf_ind1 = find(lp(:,1)<level_timing(split_level_ind,2),1,'last');
        if and(tm_ind1-ts_ind1>10, tf_ind1-tm_ind1>10)
            [~,p,~] = kstest2(lp(ts_ind1:tm_ind1,2),lp(tm_ind1:tf_ind1,2));
        else
            p = 1; % skip this level, it is too short
        end
        %display(['inds = [' num2str(ts_ind1) ', ' num2str(tm_ind1) ', ' num2str(tf_ind1) '],   p = ' num2str(p)])

        if p<pValue

            % split off levels
            if split_level_ind==1
                level_timing = [level_timing(1,:); level_timing];
                level_currents = [level_currents(1); level_currents];
            elseif split_level_ind==numel(level_currents)
                level_timing = [level_timing; level_timing(end,:)];
                level_currents = [level_currents; level_currents(end)];
            else
                level_timing = [level_timing(1:split_level_ind,:); level_timing(split_level_ind,:); level_timing((split_level_ind+1):end,:)];
                level_currents = [level_currents(1:split_level_ind); level_currents(split_level_ind); level_currents((split_level_ind+1):end)];
            end
            level_timing(split_level_ind,2) = split_time;
            level_timing((split_level_ind+1),1) = split_time;

            % update current of new level 1
            ts_ind = find(time>level_timing(split_level_ind,1),1,'first');
            tf_ind = find(time<level_timing(split_level_ind,2),1,'last');
            level_currents(split_level_ind) = mean(current(ts_ind:tf_ind));

            % update current of new level 2
            ts_ind2 = find(time>level_timing(split_level_ind+1,1),1,'first');
            tf_ind2 = find(time<level_timing(split_level_ind+1,2),1,'last');
            level_currents(split_level_ind+1) = mean(current(ts_ind2:tf_ind2));

            discrete_current = ones(numel(time),1)*level_currents(1);
            for j = 2:numel(level_currents)
                ind = find(time>level_timing(j,1));
                discrete_current(ind:end,1) = level_currents(j);
            end

        end

        % calculate the square distance of the discretized data from
        % the real data (per point)
        num = num+1;
        metric(num) = sum((discrete_current - current).^2)/numel(time);
        %f = sort((discrete_current - current).^2,'descend');
        %metric(num) = mean(f(1:100));
        fprintf('.')

    end
    
    figure(5)
    plot(1:numel(metric),metric,'o-')
    ylabel('metric')
    xlabel('number of levels')
    
    fprintf('\n')
    
    %% package results
    
    % get standard deviations of levels
    level_stds = nan(numel(level_currents),1);
    for i = 1:numel(level_currents)
        
        indices = level_timing(i,1)/sigdata.si : level_timing(i,2)/sigdata.si;
        selection = unique(sort(ceil(numel(indices)*rand(1,1e5)))); % random indices, maximum of 1e5
        level_data = sigdata.get(indices(selection),2)*1000;
        level_stds(i) = std(level_data);
        clear level_data
        fprintf('.')
        
    end
    
    fprintf('\n')
    
    discreteData.time = time;
    discreteData.current = current;
    discreteData.level_means = level_currents;
    discreteData.level_stds = level_stds;
    discreteData.level_timing = level_timing;
    
end

