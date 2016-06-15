function events = cusum_event_find(sigdata,tr,chan,threshold)
% cusum_event_find uses the cusum alrogithm to find "events" departing from
% the baseline
% Stephen Fleming 4/20/16
    %%
    % find mean baseline
    down = util.downsample_pointwise(sigdata,chan,tr,10000);
    m = mode(down);
    delta = std(down(down>m-0.01 & down<m+0.01));
    mean_baseline = mean(down(down<m+delta & down>m-delta));
%     display(mean_baseline)
    stdev = std(down(down<m+3*delta & down>m-3*delta));
    
    % go through data in chunks of 2^19 data points
    chunk = 2^19;
    ind = tr(1)/sigdata.si;
    end_ind = tr(2)/sigdata.si;
    display(' ')
    display('Cusum event finding')
    events = cell(0);
    while ind < end_ind
    
        % grab data
        ind2 = min(end_ind,ind+chunk-1);
        d = sigdata.get(ind:ind2,chan);
        
        % flip and zero signal
        d = (d-mean_baseline)*-1;
        
        % do cusum algorithm
        c(1) = 0;
        locs = [];
        for i = 2:numel(d)
            c(i) = max(0, c(i-1) + d(i) - stdev - max(0,c(i-1)-threshold*stdev)/10);
            if c(i-1)<threshold*stdev && c(i)>=threshold*stdev
                locs(end+1) = i;
            end
        end
        
        % identify peaks and then find starts and ends and means
        for i = 1:numel(locs)
            i1 = sigdata.findPrev(@(x) x(:,2)>mean_baseline-stdev, ind+locs(i));
            i2 = sigdata.findNext(@(x) x(:,2)>mean_baseline-stdev, ind+locs(i));
            events{end+1}.start_ind = i1;
            events{end}.end_ind = i2;
            events{end}.start_time = i1*sigdata.si;
            events{end}.end_time = i2*sigdata.si;
            events{end}.mean = mean(sigdata.get(i1:i2,chan));
            events{end}.open = mean_baseline;
            events{end}.min = min(sigdata.get(i1:i2,chan));
        end
        
        fprintf('.');
        ind = ind2+1;
        
    end
    display(' ')
    
%     figure()
%     view = sigdata.getViewData(tr);
%     plot(view(:,1),view(:,2))
%     hold on
%     for i = 1:numel(events)
%         plot(events{i}.start_time,mean_baseline-2*stdev,'g*')
%         plot(events{i}.end_time,mean_baseline-2*stdev,'r*')
%     end
    
end