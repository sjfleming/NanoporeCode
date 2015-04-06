function filtered = downsample2(sigdata,channel,trange,frequency)
%DOWNSAMPLE does a median downsampling to 'frequency' Hertz

% Stephen Fleming, September 23, 2014

    % number of data points
    n = round((trange(2)-trange(1))/sigdata.si);
    start = round(trange(1)/sigdata.si)+1; % index
    ending = start + n; % index

    % downsample data in chunks of 500000
    filtered = [];
    numpts = 500000;
    rep = round((1/sigdata.si)/frequency); % number of original points per downsampled point
    chunks = floor(n/numpts); % number of full chunks
    if chunks ~= 0
        for i = 1:chunks % do chunks of numpts points
            fulldata = sigdata.get(start+(i-1)*numpts:start+i*numpts-1,channel); % get chunk
            filtered = [filtered accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@median)'];
            clear fulldata
        end
    end
    if mod(n,numpts)~=0
        fulldata = sigdata.get(start+chunks*numpts:ending,channel); % the last bit that's not a full chunk
        filtered = [filtered accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@median)'];
    end

end