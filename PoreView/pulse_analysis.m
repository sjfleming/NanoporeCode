function [pulses, candidates, distances] = pulse_analysis(sigdata, trange, discreteData, tolerance)
% PULSE_ANALYSIS analyzes the timing of pulses and the proximity to level
% changes, if that data is available as the 'discreteData' structure in the
% workspace. 'tolerance' should be given as [leadtime, lagtime], for
% example, [0.004 -0.001] would look for pulses up to 4ms leading a change
% and up to 1ms lagging a level change.

% Stephen Fleming
% 9/25/14, updated 5/4/15

    pulses=[]; candidates=[]; distances=[];
    if sigdata.nsigs > 2 % we have pulsing channel
        
%         display('Analyzing pulses.')
% 
%         % number of data points
%         n = round((trange(2)-trange(1))/sigdata.si);
%         start = round(trange(1)/sigdata.si)+1; % index
%         ending = start + n; % index
% 
%         % get data in chunks, no downsampling...
%         pulseInd = [];
%         numpts = 2e6;
%         chunks = ceil((ending-start)/numpts);
%         for chunk=0:1:(chunks-1)
%             s = start + chunk*numpts;
%             pulseData = sigdata.get(s:min(s+numpts-1,ending),4);
%             if sum(pulseData)/numpts > 0.012 % skip the search if there are no pulses
%                 ind = find(pulseData<0.1,1,'first'); % if we start on a TTL high, go past it
%                 j = ind; % j goes from 1 to numpts
%                 while and(j+s < ending, ~isempty(ind))
%                     ind = find(pulseData(j+1:end)>0.5,1,'first'); % leading edge
%                     j = j+ind;
%                     pulseInd = [pulseInd j+s]; % this index signifies pulse timing
%                     j = j+round(0.001/sigdata.si)+50; % gets us past trailing edge, since TTL high lasts ~1ms
%                     display(j+s)
%                     display(ending)
%                 end
%             end
%         end
%         pulses = pulseInd*sigdata.si; % convert to time

        display('Finding pulse timings...')
        % find possible pulses quickly using a derivative
        pds = 1000;
        [pulsedata, ~] = downsample(sigdata,4,trange,pds*10,pds);
        difference = diff(pulsedata);
        [~,candidates] = findpeaks(difference.*(difference>0),'MinPeakHeight',2,'MinPeakDist',50);
        candidates = trange(1) + candidates/pds; % convert from index to time
        % refine pulse timing
        pulses = [];
        for i = 1:numel(candidates)
            ind = (candidates(i)-10e-3)/sigdata.si;
            pulses(end+1) = sigdata.findNext(@(x) x(:,4)>1, ind);
            pulses(end) = pulses(end) * sigdata.si; % convert from index to time
        end
        
        display('Finding pulses that coincide with level changes...')
        
        % find pulses that coincide with level changes
        if ~isempty(discreteData)
            delta = [];
            ind = [];
            for i=1:numel(discreteData.levelTiming(:,1))
                diffs = discreteData.levelTiming(i,1)-pulses;
                if sum(diffs>tolerance(2))~=0
                    [delta(end+1),ind(end+1)] = min(diffs(diffs>tolerance(2)));
                end
            end
            % return the ones that are within the tolerances
            criteria = and((delta>tolerance(2)),(delta<tolerance(1)));
            candidates = pulses(ind(criteria));
            distances = delta(criteria);
        end
        
        display('Done.')

    else
        display('No pulse data recorded!')
    end

end