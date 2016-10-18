classdef analysis < handle
    
    properties
        
        sigdata = [];
        voltage_view = NaN;
        current_view = NaN;
        tr = [-1,-1];
        
    end
    
    methods
        
        % constructor
        function obj = analysis(sigdata)
            % initialize the analysis object by linking it to one 
            % SignalData object
            
            obj.sigdata = sigdata;
            
        end
        
        function [cond, cond_std] = getOpenPoreConductance(obj, varargin)
            % return open pore conductance
            % based on a histogram
            % optional inputs: trange, voltage, mincond, maxcond
            
            % input handling
            in = obj.parseOptionalInputs(varargin{:});
            
            % grab conductance trace
            [voltage_raw, current_raw] = obj.getViewData(in.trange);
            voltagelogic = obj.findSpecifiedVoltageRegions(in.trange, in.voltage);
            
            % do histogram
            current = medfilt1(abs(current_raw(voltagelogic)),10);
            voltage = medfilt1(abs(voltage_raw(voltagelogic)),10);
            conductance = current ./ voltage;
            lowCond = max(in.mincond,min(conductance));
            highCond = min(in.maxcond,max(conductance));
            del = 0.002; % in nS
            xcond = lowCond-0.1:del:highCond+0.1;
            hcond = hist(conductance(voltage>5 & conductance>lowCond & conductance<highCond),xcond);
            hcond = hcond/sum(hcond);
            
            % find largest-conductance local maximum in histogram
            [~,indMax] = find(hcond>0.05,1,'last');
            [~,indMin] = find(hcond(1:indMax)<0.05,1,'last');
            [~,ind] = max(hcond(indMin:indMax));
            cond = xcond(ind+indMin-1);
            condfull = current_raw(voltagelogic)./voltage;
            cond_std = std(condfull(condfull>xcond(indMin-10*ind) & condfull<xcond(indMax+10*ind)));
        end
        
        function regions = findEventRegions(obj, varargin)
            % find events in data specified by user parameters
            % optional inputs: trange, voltage, threshhold, eventstart, mincond, maxcond
            
            % inputs
            in = obj.parseOptionalInputs(varargin{:});
            
            % get the open pore conductance
            [g_m, g_s] = obj.getOpenPoreConductance('mincond', in.mincond, 'maxcond', in.maxcond, 'voltage', in.voltage);
            
            % coarse event finding using a threshhold
            [voltage, current] = obj.getViewData(in.trange);
            dt = diff(in.trange)/numel(current);
            voltagelogic = obj.findSpecifiedVoltageRegions(in.trange, in.voltage);
            voltage = medfilt1(voltage, 10); % limit our analysis to sections with specified voltage(s)
            voltage(~voltagelogic) = NaN;
            current(~voltagelogic) = NaN;
            conductance = current./voltage;
            clear current;
            v_with_regions_deleted = voltage(voltagelogic);
            clear voltagelogic;
            V = mode(round(v_with_regions_deleted(v_with_regions_deleted>nanmax(v_with_regions_deleted)/2))); % capture voltage assumed to be most prevalent overall high voltage value
            clear v_with_regions_deleted;
            if strcmp(in.eventstart, 'currentdrop')
                lowcond = conductance < g_m * in.threshhold; % regions of conductance below threshhold
                dlogic = diff([1; lowcond; 0]);
                startcondition = @(x) x(:,2)*1000./x(:,3) < g_m * in.threshhold;
                startcondition_inv = @(x) x(:,2)*1000./x(:,3) > g_m * in.threshhold;
                endcondition = @(x,evtcond) x(:,2)*1000./x(:,3) < g_m * in.threshhold;
            elseif strcmp(in.eventstart, 'voltagedrop')
                lowvolt = voltage < V * in.threshhold; % regions of voltage below threshhold
                dlogic = diff([1; lowvolt; 0]);
                startcondition = @(x) x(:,3) < V * in.threshhold;
                startcondition_inv = @(x) x(:,3) > V * in.threshhold;
                endcondition = @(x,evtcond) x(:,2)*1000./x(:,3) < g_m * in.threshhold;
            end
            clear voltage;
            [~,possibleStartInds] = findpeaks(double(dlogic > 0),'minpeakheight',0.5,'minpeakdist',5);
            % check to make sure current starts at open pore level
            startsHigh = arrayfun(@(x) nanmax(conductance(x-5:x)) > g_m * in.threshhold, possibleStartInds);
            possibleStartInds = possibleStartInds(startsHigh);
            % one end for each start
            possibleEndInds = arrayfun(@(x) find(or(conductance(x+2:end) > g_m * in.threshhold, isnan(conductance(x+2:end))), 1, 'first'), possibleStartInds) + possibleStartInds + 1;
            % to get rid of the off-by-one errors
            possibleStartInds = possibleStartInds - 1;
            
            % exact start and end search
            start_inds = -1*ones(numel(possibleStartInds),1);
            end_inds = start_inds;
            pad = -2;
            for i = 1:numel(possibleStartInds) % go through all candidates
                % use a double-finding scheme to be as robust as possible
                % find start
                temp_start = obj.sigdata.findPrev(@(x) startcondition_inv(x), ...
                    (in.trange(1)/dt + min(possibleStartInds(i)+pad, possibleEndInds(i))) * dt/obj.sigdata.si);
                start_inds(i) = obj.sigdata.findNext(@(x) startcondition(x), temp_start-1);
                % get mean event conductance
                evtconductancearray = conductance(possibleStartInds:possibleEndInds);
                evtconductance = mean( evtconductancearray(evtconductancearray < mean([g_m, min(evtconductancearray)])) );
                % find end by next cross of threshhold
                if possibleEndInds(i)-possibleStartInds(i) < pad*dt/obj.sigdata.si
                    % end is so close we can search from the start
                    end_inds_thresh = obj.sigdata.findNext(@(x) or(x(:,2)*1000./x(:,3) > g_m * in.threshhold, x(:,3)<=0), start_inds(i)+pad);
                    end_inds_thresh = obj.sigdata.findPrev(@(x) or(endcondition(x, evtconductance), x(:,3)<=0), end_inds_thresh+1);
                else
                    % end is far enough that we should work from our
                    % best guess of the end itself
                    end_inds_thresh = obj.sigdata.findNext(@(x) or(x(:,2)*1000./x(:,3) > g_m * in.threshhold, x(:,3)<=0), ...
                        (in.trange(1)/dt + max(possibleEndInds(i)-pad, possibleStartInds(i))) * dt/obj.sigdata.si+1);
                    end_inds_thresh = obj.sigdata.findPrev(@(x) or(endcondition(x, evtconductance), x(:,3)<=0), end_inds_thresh+1);
                end
                % end should actually be when current starts to return
                % to open pore
                ending_bit = obj.sigdata.get(max(start_inds(i),end_inds_thresh-20):max(start_inds(i),end_inds_thresh-5),2)*1000;
                end_cond = mean(ending_bit);
                end_cond_std = std(ending_bit);
                end_inds(i) = obj.sigdata.findPrev(@(x) x(:,2)*1000./x(:,3) < end_cond + 0.2*end_cond_std, end_inds_thresh);
                % make sure we get some ending, if that technique
                % didn't work
                if isempty(end_inds(i))
                    end_inds(i) = end_inds_thresh;
                    display('problem identifying exact event end')
                end
                % make sure this doesn't overlap previous event
                if i>1
                    if start_inds(i)<end_inds(i-1)
                        start_inds(i) = NaN;
                        end_inds(i) = NaN;
                    end
                end
            end
            start_inds = start_inds(~isnan(start_inds));
            end_inds = end_inds(~isnan(end_inds));
            
            % give the indices of the regions as output
            regions = [start_inds, end_inds]-1; % fix off-by-one from 'find'
            %regions = [possibleStartInds*dt/obj.sigdata.si, possibleEndInds*dt/obj.sigdata.si]+obj.tr(1)/obj.sigdata.si;
        end
        
        function showEventsInPoreView(obj, pv, r, how)
            % necessary inputs: PoreView object, result of
            % findEventRegions, how ('all' or 'one')
            % plots the events, one-by-one, in PoreView
            pv.clearAxes();
            if strcmp(how,'all')
                for i = 1:size(r,1)
                    cond(i,1) = pv.data.get(r(i,1),2)*1000/pv.data.get(r(i,1),3);
                    cond(i,2) = pv.data.get(r(i,2),2)*1000/pv.data.get(r(i,2),3);
                end
                plot(r(:,1)*pv.data.si,cond(:,1),'go')
                plot(r(:,2)*pv.data.si,cond(:,2),'rx')
            elseif strcmp(how,'one')
                for i = 1:size(r,1)
                    pv.setView(max(0,sort(r(i,:)+[-100, 100])).*pv.data.si);
                    plot(r(i,1)*pv.data.si,pv.data.get(r(i,1),2)*1000/pv.data.get(r(i,1),3),'go')
                    plot(r(i,2)*pv.data.si,pv.data.get(r(i,2),2)*1000/pv.data.get(r(i,2),3),'rx')
                    pause();
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function in = parseOptionalInputs(obj, varargin)
            % parse all inputs so all methods can use them easily
            p = inputParser;
            
            % defaults and checks
            defaultFilterFreq = 1000;
            defaultSampleFreq = 5000;
            checkFilterFreq = @(x) all([isnumeric(x), numel(x)==1, x>=10, x<=20000]);
            
            defaultTimeRange = [obj.sigdata.tstart, obj.sigdata.tend]; % start and end of file
            checkTimeRange = @(x) all([all(isnumeric(x)), numel(x)==2, x(1)<x(2), x(1)>=obj.sigdata.tstart, x(2)<=obj.sigdata.tend]);
            
            checkPosNum = @(x) all([all(isnumeric(x)), all(x>=0)]);
            
            checkEventStart = @(x) any([strcmp(x, 'voltagedrop'), strcmp(x, 'currentdrop')]);
            
            % set up the inputs
            addOptional(p, 'filter', defaultFilterFreq, checkFilterFreq); % filter frequency
            addOptional(p, 'sample', defaultSampleFreq, checkFilterFreq); % (down-) sampling frequency
            addOptional(p, 'trange', defaultTimeRange, checkTimeRange); % time range of interest
            addOptional(p, 'mincond', 1.4, checkPosNum); % min open pore conductance
            addOptional(p, 'maxcond', 3, checkPosNum); % max open pore conductance
            addOptional(p, 'minduration', 1e-4, checkPosNum); % min event duration
            addOptional(p, 'threshhold', 0.95, checkPosNum); % fraction of open pore event threshhold
            addOptional(p, 'voltage', [], checkPosNum); % voltage(s) of interest
            addOptional(p, 'eventstart', 'currentdrop', checkEventStart); % what defines start of event
            
            % parse
            parse(p,varargin{:});
            in = p.Results;
        end
        
        function [voltage, current] = getViewData(obj, trange)
            % get the downsampled view data
            % only go back to file if we haven't done this already
            
            if any([numel(obj.voltage_view)==1, numel(obj.current_view)==1, all(trange ~= obj.tr)])
                raw = obj.sigdata.getViewData(trange); % grab downsampled data
                obj.voltage_view = raw(:,3);
                obj.current_view = raw(:,2)*1000;
                obj.tr = trange;
            end
            
            voltage = obj.voltage_view;
            current = obj.current_view;
            
        end
        
        function voltagelogic = findSpecifiedVoltageRegions(obj, trange, V)
            % find regions of data with specified voltage(s)
            
            % grab downsampled voltage data
            [voltage, ~] = obj.getViewData(trange);
            
            if isempty(V)
                % no voltage(s) specified
                voltagelogic = true(size(voltage));
            else
                % window of 2mV around any specified voltage of interest
                voltagelogic = any(cell2mat(arrayfun(@(x) voltage>x-2 & voltage<x+2, V, 'uniformoutput', false)), 2);
            end
            
        end
        
    end
    
end