function [filtdata] = filt_severe_median_downsample(d, n, r)
%FILT_SEVERE_MEDIAN_DOWNSAMPLE is a filter for the data which applies
%   n iterations of median filtering between downsampling at 10x or less
%   until data is downsampled by a factor of r.
%   Stephen Fleming, June 23, 2015
    
    time = d(:,1)';
    current = d(:,2)';
    
    % figure out rounds of downsampling
    while r(end)>10
        r(end+1) = r(end)/10;
        r(end-1) = 10;
    end
    display(num2str(r))
    
    for i = 1:numel(r)
        
        display(['Round ' num2str(i) ' of filtering...'])
        
        % decimate by factor of r(i)
        current = filt_decimate(current,r(i));
        time = filt_decimate(time,r(i));
        
        % median filter
        med = medfilt1(current,51,1e5);
        for j = 1:(n-1)
            med = medfilt1(med,51,1e5);
        end
        
        current = med;
        
    end
    
    % median filter
    med = medfilt1(current,111,1e5);
    for j = 1:(n-1)
        med = medfilt1(med,111,1e5);
    end
    
    filtdata = [time med];

end
