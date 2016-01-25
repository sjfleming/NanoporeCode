function [filtdata] = find_discrete_levels(sigdata, trange, cutoff)
%FIND_DISCRETE_LEVELS(trange, cutoff) uses
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.  Starts from an assumption of one level and adds
%   levels based on the size of the derivative of severely median-filtered
%   data.
%   cutoff = 0.1 works pretty well usually
%   Stephen Fleming, August 3, 2015
    
    %% divide data up into sections of 100 seconds maximum
    for i = 1:ceil(diff(trange)/100)
        
        data = sigdata.getByTime(tr);
        raw = sigdata.getViewData(tr);
        time = data(:,1);
        current = data(:,2)*1000;
        d = [time, current];
        clear data current time
        
        %% filter heavily and downsample
        
        medianRepeats = 5;
        targetFs = 200; % Hz, final frequency after downsampling
        lp = filt_lpb(d,4,1000);
        [filtered] = filt_severe_median_downsample(lp,medianRepeats,(1/pv.data.si)/targetFs);
        time = filtered(:,1);
        current = filtered(:,2);
        clear d filtered
        
        %% find level changes
        
        % take derivative
        smoothing = 15;
        deriv = smooth_derivative(current, smoothing);
        
        % start up the arrays of currents and times
        level_currents = mean(current);
        level_timing = [time(1) time(end)];
        pValue = 1e-5;
        
        % find peaks in the derivative and sort them
        [deriv_peaks, deriv_peak_locs] = findpeaks(abs(deriv),'MinPeakDist',2*smoothing,'MinPeakHeight',0.25*std(deriv));
        [~, sort_ind] = sort(deriv_peaks,'descend');
        num = 1;
        metric = 100;
        
        while metric(end)>cutoff
        %for i = 1:numel(deriv_peaks)
            
            % find where level would split
            split_time = time(deriv_peak_locs(sort_ind(num)));
            split_level_ind = find(split_time>level_timing(:,1),1,'last');
            
            % check to make sure levels are statistically different
            ts_ind1 = find(time>level_timing(split_level_ind,1),1,'first');
            tf_ind1 = find(time<level_timing(split_level_ind,2),1,'last');
            [~,p,~] = kstest2(current(ts_ind1:deriv_peak_locs(sort_ind(num))),current(deriv_peak_locs(sort_ind(num)):tf_ind1));
            
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
            
        end
        
        %%
        
        optimumLevels = num;
        
    end
    
end

