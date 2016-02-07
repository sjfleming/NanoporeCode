function [filtered, unfiltered] = downsample1(sigdata,channel,trange,filter,frequency)
%DOWNSAMPLE does a median downsampling to 'frequency' Hertz, after low-pass
%   filtering at 'filter' Hertz.

% Stephen Fleming, September 23, 2014

    % number of data points
    n = round((trange(2)-trange(1))/sigdata.si);
    start = round(trange(1)/sigdata.si)+1; % index
    ending = start + n; % index

    % downsample data in chunks of 500000
    filtered = [];
    unfiltered = [];
    numpts = 500000;
    rep = round((1/sigdata.si)/frequency); % number of original points per downsampled point
    chunks = floor(n/numpts); % number of full chunks
    if chunks ~= 0
        for i = 1:chunks % do chunks of numpts points
            fulldata = sigdata.get(start+(i-1)*numpts:start+i*numpts-1,channel); % get chunk
            if channel==2
                fulldata = fulldata*1000; % put in pA
            end
            % filter the raw data first
            wn = filter/((1/sigdata.si)/2); % filter at 1kHz
            %[b,a] = butter(4, wn, 'low');
            [b,a] = maxflat(10,2,wn);
            filt = filtfilt(b,a,fulldata);
            filtered = [filtered accumarray(1+floor((1:numel(filt))/rep)',filt',[],@median)'];
            unfiltered = [unfiltered accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@median)'];
            clear fulldata
            clear filt
        end
    end
    if mod(n,numpts)~=0
        fulldata = sigdata.get(start+chunks*numpts:ending,channel); % the last bit that's not a full chunk
        if channel==2
            fulldata = fulldata*1000; % put in pA
        end
        % filter the raw data first
        wn = filter/((1/sigdata.si)/2); % filter at 'filter' Hz
        [b,a] = butter(4, wn, 'low');
        filt = filtfilt(b,a,fulldata);
        filtered = [filtered accumarray(1+floor((1:numel(filt))/rep)',filt',[],@median)'];
        unfiltered = [unfiltered accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@median)'];
    end

end