function transition = next_transition(sigdata,channel,trange,filter,minHeight)
%NEXT_TRANSITION loops through data in one second chunks and finds the
%   time of the first transition to a new level
%   Stephen Fleming
    
    % loop until point is found or end is reached
    transition = [];
    found = false;
    t = trange(1);
    step = 2;
    while and( found==false, t<trange(2) )
        
        % generate candidate transition points quickly
        frequency = 10*filter;
        [filtered,~] = downsample(sigdata,channel,[t,min(t+step,trange(2))],filter,frequency);
        filtered = medfilt1(filtered,11,1e5);
        diff = filtered(4:end)-filtered(1:end-3);
%         raw = sigdata.getViewData([t,min(t+2,trange(2))]);
%         rawSamplingInterval = raw(2,1)-raw(1,1);
%         raw = medfilt1(raw(:,channel),51,1e5);
%         diff = raw(4:end)-raw(1:end-3);
        if(minHeight>0)
            diff = diff.*(diff>0);
        else % flip so peaks we want are positive
            diff = -1*diff.*(diff<0);
            minHeight = -1*minHeight;
        end
%         figure(2)
%         plot(linspace(t,min(t+2,trange(2)),numel(diff)),diff)
%         pause()
        [~,candidates] = findpeaks(diff,'MinPeakHeight',minHeight,'MinPeakDistance',50);
       candidates = candidates*((1/sigdata.si)/frequency); % translate candidate indices to original sampling frequency
%         candidates = candidates*(rawSamplingInterval/sigdata.si); % translate candidate indices to original sampling frequency
        
        % check candidates
        tol = 1e-15;
        startInd = t*(1/sigdata.si);
        if ~isempty(candidates)
            for j = 1:numel(candidates)
                % find exact location of transition
                num = 41;
                pValues = ones(1,num);
                data = sigdata.get(startInd+candidates(j)-100:startInd+candidates(j)+100 ,channel);
                for i = -1*(num-1):(num-1)
                    %before = sigdata.get( startInd+j-100:startInd+j+i-1 ,channel);
                    %after = sigdata.get( startInd+j+i+1:startInd+j+100 ,channel);
                    before = data(1:100+i);
                    after = data(102+i:end);
                    [~,pValues(num+i)] = kstest2(before,after);
                end
                % break on first one we find
                [lowestP,ind] = min(pValues);
                if lowestP<tol
                    transition = t+(candidates(j)+(ind-num))*sigdata.si; % time in seconds
                    found = true;
                    break;
                end
            end
        end
        t = t+step;
        
    end
    
end