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
    
    fprintf('\b')
    
    % iteratively downsample and median filter (multiple times)
    for i = 1:numel(r)
        
        % decimate by factor of r(i)
        current = filt_decimate(current,round(r(i)));
        %current = decimate(current,round(r(i)))';
        fprintf('.')
        
        % median filter
        med = medfilt1(current,51,1e5);
        for j = 1:(n-1)
            med = medfilt1(med,51,1e5);
            fprintf('.')
        end
        
        current = med;
        
    end
    
%     % median filter
%     med = medfilt1(current,111,1e5);
%     for j = 1:(n-1)
%         med = medfilt1(med,111,1e5);
%         fprintf('.')
%     end
    
    time = linspace(time(1),time(end),numel(med))';
    filtdata = [time med];
    fprintf('\n')

end
