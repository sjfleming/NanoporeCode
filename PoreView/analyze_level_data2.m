function [discreteData, V] = analyze_level_data2(sigdata, trange, chan, fs)
%ANALYZE_LEVEL_DATA loops through data between cursors (or all) and applies
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.
%   Stephen Fleming, March 5, 2015
    
    % Downsample to fs
    display('Finding discrete levels:')
    currentData = downsample2(sigdata,chan,trange,fs)*1000; % in pA
    time = linspace(trange(1),trange(2),numel(currentData));
    
    % Find possible transitions
    display('Identifying level transitions...')
    n = 10;
    for i = 2:numel(currentData)-1 % smoothing
        before = mean(currentData(max(i-n,1):i-1));
        after = mean(currentData(i+1:min(i+n,numel(currentData))));
        diff(i) = after-before;
    end
    [~,index] = findpeaks(abs(diff).*(abs(diff)-0.3*std(diff)>0),'minpeakdistance',fs*0.01);
    
    % Eliminate unlikely transitions using Kolmogorov-Smirnov test
    display('Refining...')
    inds = [1 index numel(currentData)];
    worsts = [numel(inds) 1];
    pValues = [inds' zeros(numel(inds),1)];
    for i = 2:numel(inds)-1
        [~,pValues(i,2)] = kstest2(currentData(inds(i-1):inds(i)),currentData(inds(i):inds(i+1)));
    end
    maxlevels = 5000; % set a hard cutoff for maximum number of levels that can be found
    % pCutoff = 1e-80; % good for data with extreme 1/f noise
    % pCutoff = 1e-25; % good for low-noise data
    pCutoff = 1e-50; % good for med-noise data
    while or(and(worsts(end,2) > pCutoff, size(pValues,1) > 2), size(pValues,1) > maxlevels) % iterate until the largest p value falls below cutoff
        [pBad,indsind] = max(pValues(:,2)); % find worst p-value
        worsts(end+1,:) = [numel(inds) pBad]; % keep track of it
        inds(indsind) = [];
        pValues(indsind,:) = [];
        if indsind == 2
            first = indsind;
        elseif indsind ~=1
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

