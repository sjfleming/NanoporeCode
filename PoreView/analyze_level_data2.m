function [discreteData, V] = analyze_level_data2(sigdata, trange, chan, fs, t)
%ANALYZE_LEVEL_DATA loops through data between cursors (or all) and applies
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.
%   Stephen Fleming, March 5, 2015
    
    % Downsample to fs
    display(['Plotting time interval [' num2str(trange) ']'])
    currentData = downsample2(sigdata,chan,trange,fs)*1000; % in pA
    time = linspace(trange(1),trange(2),numel(currentData));
    
    %% make pretty plot

    % Get the reduced version of the unfiltered data as well
    raw = sigdata.getViewData(trange);
    
    % Plot!
    figure(2)
    clf(2)
    h = gcf;
    plot(raw(:,1),raw(:,2)*1000,'Color',[0.85,0.85,0.85])
    hold on
    if sigdata.nsigs>2 % we have pulses recorded
        [pulses,~,~] = pulse_analysis(sigdata,trange,[],[]); % get pulse timings
        line(repmat(pulses',1,2)',repmat(get(gca,'ylim'),numel(pulses),1)','Color',[0.9 1 0.9],'LineStyle','-') % vertical lines
    end
    plot(time,currentData,'Color',[0,0.4470,0.7410])
    xlabel('Time (s)','FontSize',28)
    ylabel('Current (pA)','FontSize',28)
    name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];
    
    title(t,'FontSize',24)
    annotation('textbox', [0.75 0.87 0 0], 'String', name, 'FontSize', 20);
    xlim([trange(1), trange(2)])
    ylim([min(0,min(currentData)-20),max(currentData)+20])
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
    set(h,'Position',[100 500 1000 400]) % size the figure
    set(gca,'FontSize',24)
    
    %% Find possible transitions
    display('Identifying level transitions...')
    n = 10;
%     for i = 2:numel(currentData)-1 % smoothing
%         before = mean(currentData(max(i-n,1):i-1));
%         after = mean(currentData(i+1:min(i+n,numel(currentData))));
%         diff(i) = after-before;
%     end
%     bigData = nan(n,numel(currentData)+n-1);
%     for i = 1:n
%         bigData(i,i:(i+numel(currentData)-1)) = currentData;
%     end
%     bigData2 = nanmean(bigData,1);
    difference = diff(currentData);
    [~,index] = findpeaks(abs(difference).*(abs(difference)-0.3*std(difference)>0),'minpeakdistance',fs*0.01);
    
    % Eliminate unlikely transitions using Kolmogorov-Smirnov test
    display('Refining...')
    inds = [1 index numel(currentData)];
    worsts = [numel(inds) 1];
    pValues = [inds' zeros(numel(inds),1)];
    for i = 2:numel(inds)-1
        [~,pValues(i,2)] = kstest2(currentData(inds(i-1):inds(i)),currentData(inds(i):inds(i+1)));
    end
    maxlevels = 3000; % set a hard cutoff for maximum number of levels that can be found
    %pCutoff = 1e-20; % good for low-noise data
    pCutoff = 1e-80; % good for med-noise data
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
        levels_std(i) = std(currentData(inds(i):inds(i+1)));
        current(inds(i):inds(i+1)) = levels(i)*ones(1,inds(i+1)-inds(i)+1);
        levelTiming(i,:) = [time(inds(i)) time(inds(i+1))];
    end
    
    % Get applied voltage
    V = round(sigdata.get(trange(1)/sigdata.si,3)/10)*10; % applied voltage in mV
    
    % Set variables to pass back
    discreteData.current = current;
    discreteData.time = time;
    discreteData.levels = levels';
    discreteData.levels_std = levels_std';
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

