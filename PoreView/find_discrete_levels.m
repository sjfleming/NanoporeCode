function discreteData = find_discrete_levels(sigdata, channel, trange, finalFrequency, pValue)
% discreteData = find_discrete_levels(sigdata, trange, finalFrequency, medianRepeats, pValue)
%   downsamples and median filters heavily in chunks of 100 sec max, and
%   then runs Laszlo's level-finding algorithm
%   Stephen Fleming, October 7, 2015
    
    %% divide data up into sections of 100 seconds maximum, filter heavily and downsample
    time = [];
    current = [];
    for i = 1:ceil(diff(trange)/100)
        
        tr = [trange(1)+(i-1)*100, min(trange(2),trange(1)+i*100)];
        display(['Working on time range [' num2str(tr(1)) ', ' num2str(tr(2)) ']'])
        data = sigdata.getByTime(tr);
        t = data(:,1);
        c = data(:,channel)*1000;
        d = [t, c];
        clear data c t
        
        [filtered] = filt_decimate(d,(1/sigdata.si)/finalFrequency);
        %filtered = filt_med(filtered,11);
        t = filtered(:,1);
        c = filtered(:,2);
        clear d filtered
        
        time = [time; t];
        current = [current; c];
        
    end
    
    %% use Laszlo's level-finding algorithm to find levels
    
    levels = laszlo_levels([time, current],pValue);
    
    %% package the data
    
    discreteData.time = time;
    discreteData.current = current;
    discreteData.level_means = cellfun(@(x) x.current_mean, levels);
    discreteData.level_medians = cellfun(@(x) x.current_median, levels);
    discreteData.level_stds = cellfun(@(x) x.current_std, levels);
    discreteData.level_timing = [cellfun(@(x) x.start_time, levels), cellfun(@(x) x.end_time, levels)];
    
end

