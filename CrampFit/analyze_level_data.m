function [discreteData, V] = analyze_level_data(sigdata, trange)
%ANALYZE_LEVEL_DATA loops through data between cursors (or all) and applies
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.
%   Stephen Fleming, September 23, 2014
    
    % Low pass filter at 1kHz and downsample to 10kHz
    display('Finding discrete levels:')
    factor = min(10,ceil((trange(2)-trange(1))/20));
    f = 1000/factor;
    fs = 10000/factor;
    display(['Filtering at ' num2str(f) 'Hz and downsampling to ' num2str(fs) 'Hz...'])
    [currentData, unfilteredCurrent] = downsample(sigdata,2,trange,f,fs);
    time = linspace(trange(1),trange(2),numel(currentData));
    
    % Median filter
    display('Median filtering heavily...')
    i = 1;
    delta = 10;
    filterData = medfilt1(currentData,51,1e5);
    while delta>0.0001
        filterData2 = medfilt1(filterData,51,1e5); % window is ~2ms: longer than one level, shorter than two
        delta = sum(abs(filterData2-filterData));
        i = i+1;
        filterData = filterData2;
    end
    
    % Find possible transitions
    display('Identifying level transitions...')
    n = 10;
    for i = 2:numel(filterData)-1 % smoothing
        before = mean(filterData(max(i-n,1):i-1));
        after = mean(filterData(i+1:min(i+n,numel(filterData))));
        diff(i) = after-before;
    end
    [~,index] = findpeaks(abs(diff).*(abs(diff)-0.3*std(diff)>0),'minpeakdistance',fs*0.01);
%     figure()
%     plot(abs(diff))
%     hold on
%     line([1,numel(diff)],1*[std(diff) std(diff)])
    
    % Eliminate unlikely transitions using Kolmogorov-Smirnov test
    display('Refining...')
    inds = [1 index numel(currentData)];
    worsts = [numel(inds) 1];
    pValues = [inds' zeros(numel(inds),1)];
    for i = 2:numel(inds)-1
        [~,pValues(i,2)] = kstest2(unfilteredCurrent(inds(i-1):inds(i)),unfilteredCurrent(inds(i):inds(i+1)));
    end
    while worsts(end,2) > 1e-15 && size(pValues,1) > 2 % iterate until the largest p value falls below cutoff
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
            %display(['Recalculating indices ' num2str(first) ' and ' num2str(second)])
            for i = first:second  % recalculate the ones before and after (if there is a before and after)
                [~,pValues(i,2)] = kstest2(unfilteredCurrent(inds(i-1):inds(i)),unfilteredCurrent(inds(i):inds(i+1)));
            end
        end
    end
    for i = 1:numel(inds)-1 % find the levels and save timing
        levels(i) = mean(currentData(inds(i):inds(i+1)));
        current(inds(i):inds(i+1)) = levels(i)*ones(1,inds(i+1)-inds(i)+1);
        levelTiming(i,:) = [time(inds(i)) time(inds(i+1))];
    end
    
    % Get applied voltage
    V = round(sigdata.get(trange(1)/sigdata.si,3)/10)*10; % applied voltage in mV
    
    % Set variables to pass back
    discreteData.current = current;
    discreteData.time = time;
    discreteData.levels = levels';
    discreteData.levelTiming = levelTiming;
    
    % Refine the timing of the level changes using full (not downsampled) data
    discreteData = refine_level_data(sigdata, discreteData);
    
    display('Done.')
    
    % add discrete version to pretty plot, if it's been made
    if ishandle(2)
        figure(2)
        hold on
        plot(discreteData.time,discreteData.current,'r','LineWidth',1)
    end
    
end

