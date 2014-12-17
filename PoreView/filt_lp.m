function [ filtdata ] = filt_lp( data, n, wp )
    %FILT_HP Filters data using a low-pass Butterworth filter
    %   Frequency wp is in Hz, n is poles.
    %   Data is passed with colums of [time, sig1, sig2, ...], and
    %       returns the same array but with the data filtered.


    % extract si. easier than passing each time
    si = data(2,1)-data(1,1);

    % convert from absolute frequency to 'normalized frequency'
    % which is 0 at 0 and 1 at Nyquist freq 1/(2*si)
    wn = 2*wp*si;
    
    if (wn > 1)
        filtdata = data;
        return
    end

    % create the filter coefficients
    [b a] = butter(n, wn, 'low');

    % copy
    filtdata = data;
    % and replace data portion for all dimensions
    for i=2:size(filtdata,2)
        filtdata(:,i) = filtfilt(b,a,data(:,i));
    end
end
