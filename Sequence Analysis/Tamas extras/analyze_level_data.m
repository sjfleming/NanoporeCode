function discreteData = analyze_level_data(sigdata, trange)
%ANALYZE_LEVEL_DATA loops through data between cursors (or all) and applies
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.
%   Stephen Fleming, September 23, 2014
    
    % Low pass filter at 1kHz and downsample to 10kHz
    display('Finding discrete levels:')
    display('Filtering and downsampling...')
    currentData = downsample(sigdata,trange,1000,10000);
    time = linspace(trange(1),trange(2),numel(currentData));
    
    % Median filter
    display('Median filtering heavily...')
    i = 1;
    delta = 10;
    while delta>0.0001
        currentData2 = medfilt1(currentData,51,1e5); % window is ~2ms: longer than one level, shorter than two
        delta = sum(abs(current2-currentData));
        i = i+1;
        currentData = currentData2;
    end
    
    % Find possible transitions
    display('Identifying level transitions...')
    n = 10;
    for i = 2:numel(currentData)-1 % smoothing
        before = mean(f(max(i-n,1):i-1));
        after = mean(f(i+1:min(i+n,numel(f))));
        diff(i) = after-before;
    end
    [~,index] = findpeaks(abs(diff).*(abs(diff)-3*std(diff)>0),'minpeakdistance',1000);
    
    % Eliminate unlikely transitions using Kolmogorov-Smirnov test
    display('Refining...')
    inds = [1 index numel(currentData)];
    worsts = [numel(inds) 1];
    pValues = [inds' zeros(numel(inds),1)];
    for i = 2:numel(inds)-1
        [~,pValues(i,2)] = kstest2(currentData(inds(i-1):inds(i)),currentData(inds(i):inds(i+1)));
    end
    while worsts(end,2) > 0.001 && size(pValues,1) > 2 % iterate until the largest p value falls below cutoff
        [pBad,indsind] = max(pValues(:,2)); % find worst p-value
        worsts(end+1,:) = [numel(inds) pBad]; % keep track of it
        inds(indsind) = [];
        pValues(indsind,:) = [];
        if indsind == 2
            first = indsind;
        else
            first = indsind-1;
        end
        if indsind == size(pValues,1)
            second = indsind-1;
        else
            second = indsind;
        end
        if first < second
            for i = first:second  % recalculate the ones before and after (if there is a before and after)
                [~,pValues(i,2)] = kstest2(currentData(inds(i-1):inds(i)),currentData(inds(i):inds(i+1)));
            end
        end
    end
    for i = 1:numel(inds)-1 % find the levels
        levels(i) = mean(currentData(inds(i):inds(i+1)));
        current(inds(i):inds(i+1)) = levels(i)*ones(1,inds(i+1)-inds(i)+1);
    end
    
    % Set variables to pass back
    discreteData = [time' current'];
%     discreteData.current = current;
%     discreteData.time = time;
%     discreteData.levels = levels;
%     display('Done.')
    
end

