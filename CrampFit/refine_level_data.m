function discreteData = refine_level_data(sigdata, discreteData)
%REFINE_LEVEL_DATA works with original data, no downsampling, and applies
%   a level-finding algorithm based on Komogorov-Smirnov testing to 
%   discretize the data, based on estimates already obtained with
%   downsampled data.
%   Stephen Fleming, October 6, 2014

display('Refining level timing...')
transitions = discreteData.levelTiming(1:end-1,2); % coarse level transition times
refinedTransitions = zeros(numel(transitions),1);

% load a 4ms window around these times and re-calculate best transition time
for i = 1:numel(transitions)
    data = sigdata.getByTime(transitions(i)-0.002,transitions(i)+0.002);
    t = data(:,1); % time in seconds
    current = data(:,2)*1000; % current in pA
    clear data
    pValues = ones(numel(current),1);
    for j = 50:numel(current)-50
        [~,pValues(j,1)] = kstest2(current(1:j-1),current(j+1:end));
    end
    [~,ind] = min(pValues);
    refinedTransitions(i,1) = t(ind);
    %display([num2str(i/numel(transitions)*100,3) '% complete'])
end

% Re-save the information about level timing in discreteData
discreteData.levelTiming(2:end,1) = refinedTransitions;
discreteData.levelTiming(1:end-1,2) = refinedTransitions;

end
