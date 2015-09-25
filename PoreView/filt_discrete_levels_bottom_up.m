function [filtdata] = filt_discrete_levels_bottom_up(d, pCutoff)
%FILT_DISCRETE_LEVELS is a filter for the data which applies
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data.
%   Stephen Fleming, June 16, 2015
    
    %% Find possible transitions
    display('Identifying level transitions...')
    data = d(:,2)';
    difference = diff(data);
    [~,index] = findpeaks(abs(difference).*(abs(difference)-0.2*std(difference)>0),'minpeakdistance',1000);
    
    % Eliminate unlikely transitions using Kolmogorov-Smirnov test
    display('Refining...')
    inds = [1 index numel(data)];
    worsts = [numel(inds) 1];
    pValues = [inds' ones(numel(inds),1)];
    for i = 2:numel(inds)-1
        [~,pValues(i,2)] = kstest2(data(inds(i-1):inds(i)),data(inds(i):inds(i+1)));
    end
    maxlevels = 3600; % set a hard cutoff for maximum number of levels that can be found
    
    while or(and(worsts(end,2) < pCutoff, size(pValues,1) > 2), size(pValues,1) > maxlevels) % iterate until the largest p value falls below cutoff
        [pBad,indsind] = min(pValues(:,2)); % find worst p-value
        display(pBad);
        worsts(end+1,:) = [numel(inds) pBad]; % keep track of it
        inds(indsind) = [];
        pValues(indsind,:) = [];
        if indsind == 2
            first = indsind;
        elseif indsind ~=1
            first = indsind-1;
        else
            first = 2;
        end
        if indsind == size(pValues,1)
            second = indsind-1;
        else
            second = indsind;
        end
        if first < second
            for i = first:second  % recalculate the ones before and after (if there is a before and after)
                [~,p,stat] = kstest2(data(inds(i-1):inds(i)),data(inds(i):inds(i+1)));
                %[~,pValues(i,2)] = ttest2(data(inds(i-1):inds(i)),data(inds(i):inds(i+1)));
                pValues(i,2) = stat;
            end
        end
    end
    
    % Find the levels and save timing
    for i = 1:numel(inds)-1
        levels(i) = mean(data(inds(i):inds(i+1)));
        transition(i) = d(inds(i),1);
        filtdata(inds(i):inds(i+1)) = levels(i)*ones(1,inds(i+1)-inds(i)+1);
    end
    
    assignin('base','level_transitions',transition);
    assignin('base','level_currents',levels);
    
    filtdata = [d(:,1), filtdata'];
    
end

