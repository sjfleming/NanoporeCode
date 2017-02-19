classdef analysis < handle
    
    properties
        
        sigdata = [];
        voltage_view = NaN;
        current_view = NaN;
        tr = [-1,-1];
        in = []; % this is the struct containing parsed inputs
        parsed = false;
        filename = ''; % where data is saved
        
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
            obj.parseOptionalInputs(varargin{:});
            
            % grab conductance trace
            [voltage_raw, current_raw] = obj.getViewData(obj.in.trange);
            voltagelogic = obj.findSpecifiedVoltageRegions(obj.in.trange, obj.in.voltage);
            
            % do histogram
            current = medfilt1(abs(current_raw(voltagelogic)),10);
            voltage = medfilt1(abs(voltage_raw(voltagelogic)),10);
            conductance = current ./ voltage;
            lowCond = max(obj.in.mincond,min(conductance));
            highCond = min(obj.in.maxcond,max(conductance));
            del = 0.01; % in nS
            xcond = lowCond-0.1:del:highCond+0.1;
            hcond = hist(conductance(voltage>5 & conductance>lowCond & conductance<highCond),xcond);
            hcond = hcond/sum(hcond);
            
            % find largest-conductance local maximum in histogram
            [~,indMax] = find(hcond>0.05,1,'last');
            [~,indMin] = find(hcond(1:indMax)<0.05,1,'last');
            [~,ind] = max(hcond(indMin:indMax));
            condfull = abs(current_raw(voltagelogic))./voltage;
            try
                cond_interest = condfull(condfull>xcond(indMin-5*ind) & condfull<xcond(min(numel(xcond),indMax+5*ind)));
                cond = mean(cond_interest);
                cond_std = std(cond_interest);
            catch ex
                pv_launch(obj.sigdata.filename);
                cond = input('Could not identify open pore conductance.  Please enter in nS: ');
                cond_std = 0.05;
            end
        end
        
        function voltages = getAppliedVoltages(obj, varargin)
            % return the main voltages applied in this file (>5% total time)
            % optional inputs: trange, voltagecheck function handle
            
            % input handling
            obj.parseOptionalInputs(varargin{:});
            
            % grab voltage trace
            d = obj.sigdata.getViewData(obj.in.trange);
            voltage = medfilt1(d(:,3),10); % median filtered, rough data
            voltage = voltage(obj.in.voltagecheck(voltage)); % apply the function 'voltagecheck' to see which time ranges have valid voltages
            
            % do histogram
            del = 1; % in mV
            xv = min(voltage)-10:del:max(voltage)+10;
            hv = hist(voltage,xv);
            hv = hv/sum(hv);
            
            % find peaks in histogram (minimum 5% of data)
            [~,inds] = findpeaks(hv,'minpeakheight',0.05,'minpeakdist',3);
            voltages = round(xv(inds));
        end
        
        function regions = findEventRegions(obj, varargin)
            % find events in data specified by user parameters
            % optional inputs: trange, voltage, threshold, eventstart, mincond, maxcond
            
            % inputs
            obj.parseOptionalInputs(varargin{:});
            
            % get the open pore conductance
            [g_m, g_s] = obj.getOpenPoreConductance(varargin{:});
            
            % coarse event finding using a threshold
            [voltage, current] = obj.getViewData(obj.in.trange);
            dt = diff(obj.in.trange)/numel(current);
            voltagelogic = obj.findSpecifiedVoltageRegions(obj.in.trange, obj.in.voltage);
            voltage = medfilt1(voltage, 10); % limit our analysis to sections with specified voltage(s)
            voltage(~voltagelogic) = NaN;
            current(~voltagelogic) = NaN;
            conductance = current./voltage;
            clear current;
            %v_with_regions_deleted = voltage(voltagelogic);
            clear voltagelogic;
            %V = mode(round(v_with_regions_deleted(v_with_regions_deleted>nanmax(v_with_regions_deleted)/2))); % capture voltage assumed to be most prevalent overall high voltage value
            clear v_with_regions_deleted;
            if strcmp(obj.in.eventstart, 'currentdrop')
                lowcond = abs(conductance) < abs(g_m) * obj.in.threshold; % regions of conductance below threshold
                dlogic = diff([1; lowcond; 0]);
                startcondition = @(x) abs(x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling)) < abs(g_m) * obj.in.threshold;
                startcondition_inv = @(x) abs(x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling)) > abs(g_m) * obj.in.threshold;
                endcondition = @(x,evtcond) abs(x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling)) < abs(g_m) * obj.in.threshold;
            elseif strcmp(obj.in.eventstart, 'voltagedrop')
                lowvolt = abs(voltage) < abs(V) * obj.in.threshold; % regions of voltage below threshold
                dlogic = diff([1; lowvolt; 0]);
                startcondition = @(x) abs(x(:,3)*obj.in.voltagescaling) < abs(V) * obj.in.threshold;
                startcondition_inv = @(x) abs(x(:,3)*obj.in.voltagescaling) > abs(V) * obj.in.threshold;
                endcondition = @(x,evtcond) abs(x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling)) < abs(g_m) * obj.in.threshold;
            end
            clear voltage;
            [~,possibleStartInds] = findpeaks(double(dlogic > 0),'minpeakheight',0.5,'minpeakdist',5);
            % check to make sure current starts at open pore level
            startsHigh = arrayfun(@(x) nanmax(conductance(max(1,x-5):x)) > g_m * obj.in.threshold, possibleStartInds);
            possibleStartInds = possibleStartInds(startsHigh);
            % one end for each start
            conductance = [conductance; nan]; % just so it won't try to go past end
            possibleEndInds = arrayfun(@(x) find(or(conductance(x+2:end) > g_m * obj.in.threshold, isnan(conductance(x+2:end))), 1, 'first'), possibleStartInds) + possibleStartInds + 1;
            % to get rid of the off-by-one errors
            possibleStartInds = possibleStartInds - 1;
            
            % trim out ones that are too short
            too_short = (possibleEndInds-possibleStartInds)*dt < obj.in.minduration;
            possibleStartInds = possibleStartInds(~too_short);
            possibleEndInds = possibleEndInds(~too_short);
            
            % exact start and end search
            start_inds = -1*ones(numel(possibleStartInds),1);
            end_inds = start_inds;
            pad = 2;
            for i = 1:numel(possibleStartInds) % go through all candidates
                % use a double-finding scheme to be as robust as possible
                % find start
                temp_start = max(2, obj.sigdata.findPrev(@(x) startcondition_inv(x), ...
                    (obj.in.trange(1)/dt + min(possibleStartInds(i)-pad, possibleEndInds(i))) * dt/obj.sigdata.si));
                start_inds(i) = obj.sigdata.findNext(@(x) startcondition(x), temp_start-1);
                % get mean event conductance
                evtconductancearray = conductance(possibleStartInds:possibleEndInds);
                evtconductance = mean( evtconductancearray(evtconductancearray < mean([g_m, min(evtconductancearray)])) );
                % find end by next cross of threshold
                start_pt = (obj.in.trange(1)/dt + max(possibleEndInds(i)-pad, possibleStartInds(i))) * dt/obj.sigdata.si+1;
                if possibleEndInds(i)-possibleStartInds(i) < pad*dt/obj.sigdata.si
                    % end is so close we can search from the start
                    start_pt = start_inds(i) + pad + round(1e-4/obj.sigdata.si);
                end
                end_inds_thresh = obj.sigdata.findNext(@(x) x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling) > g_m * obj.in.threshold | ... % open pore
                    round(x(:,3)*obj.in.voltagescaling/5)==0 | ... % voltage zero
                    x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling) < 0, ... % conductance tanked
                    start_pt);
                end_inds_thresh = obj.sigdata.findPrev(@(x) or(endcondition(x, evtconductance), round(x(:,3)*obj.in.voltagescaling/5)==0), end_inds_thresh+1);
                % end should actually be when current starts to return
                % to open pore
                ending_bit = obj.sigdata.get(max(start_inds(i),end_inds_thresh-20):max(start_inds(i),end_inds_thresh-5),2) * obj.in.currentscaling;
                V = mode(round(obj.voltage_view(max(1,possibleStartInds(i)):min(numel(obj.voltage_view),possibleEndInds(i))))); % estimate of this event's voltage
                end_cond = abs(mean(ending_bit)/V);
                %end_cond_std = std(ending_bit);
                end_inds(i) = obj.sigdata.findPrev(@(x) x(:,2)*obj.in.currentscaling./(x(:,3)*obj.in.voltagescaling) < mean([end_cond, end_cond, end_cond, g_m * obj.in.threshold]), end_inds_thresh);
                % make sure we get some ending, if that technique
                % didn't work
                if isempty(end_inds(i))
                    end_inds(i) = end_inds_thresh;
                    display('problem identifying exact event end')
                end
                % make sure this doesn't overlap previous event
                if i>1
                    if end_inds(i)<start_inds(i-1)
                        start_inds(i) = NaN;
                        end_inds(i) = NaN;
                        display('problem identifying exact event end: overlap')
                    end
                end
            end
            start_inds = start_inds(~isnan(start_inds));
            end_inds = end_inds(~isnan(end_inds));
            
            % give the indices of the regions as output
            regions = round([start_inds, end_inds]-1); % fix off-by-one from 'find'
            %regions = [possibleStartInds*dt/obj.sigdata.si, possibleEndInds*dt/obj.sigdata.si]+obj.tr(1)/obj.sigdata.si;
            
            % get rid of events that are too short
            too_short = (end_inds-start_inds)*obj.sigdata.si < obj.in.minduration;
            regions = regions(~too_short,:);
        end
        
        function showEventsInPoreView(obj, pv, events, how)
            % showEventsInPoreView(pv, events, all_or_one)
            % ex: obj.showEventsInPoreView([], event{1}, 'one');
            % necessary inputs:
            % PoreView object or []
            % result of findEventRegions or calculateEventStatistics
            % how (can be 'all' or 'one'). 'one' is for one-at-a-time, with
            % pauses in between.
            % plots the events in PoreView
            
            r = cell2mat(cellfun(@(x) x.index, events, 'uniformoutput', false));
            if isempty(pv)
                pv = pv_launch(events{1}.file);
            end
            
            pv.clearAxes();
            if strcmp(how,'all')
                for i = 1:size(r,1)
                    % use the y values for whatever signal is currently
                    % plotted
                    y(i,1) = pv.data.get(r(i,1),pv.psigs(1).sigs);
                    y(i,2) = pv.data.get(r(i,2),pv.psigs(1).sigs);
                end
                plot(pv.psigs(1).axes, r(:,1)*pv.data.si,y(:,1),'go','MarkerSize',10)
                plot(pv.psigs(1).axes, r(:,2)*pv.data.si,y(:,2),'rx','MarkerSize',10)
            elseif strcmp(how,'one')
                for i = 1:size(r,1)
                    % if PoreView's SignalData is not current, update it
                    if ~strcmp(pv.data.filename,events{i}.file)
                        pv = pv.loadFile(events{i}.file);
                    end
                    % use the y values for whatever signal is currently
                    % plotted
                    pv.setView(max(0,sort(r(i,:)+[-100, 100])).*pv.data.si);
                    plot(pv.psigs(1).axes, r(i,1)*pv.data.si, pv.data.get(r(i,1),pv.psigs(1).sigs),'go','MarkerSize',10)
                    plot(pv.psigs(1).axes, r(i,2)*pv.data.si, pv.data.get(r(i,2),pv.psigs(1).sigs),'rx','MarkerSize',10)
                    % if there are levels specified, show them
                    if isfield(events{i},'levels')
                        timing = cell2mat(cellfun(@(x) [x.start_time, x.end_time], events{i}.levels, 'uniformoutput', false));
                        means = cellfun(@(x) x.current_mean, events{i}.levels) / obj.in.currentscaling; % pA back to initial data value
                        plot(pv.psigs(1).axes, timing', (means*[1,1])', '-', 'LineWidth', 4, 'Color', 'k');
                    end
                    pause();
                end
            end
        end
        
        function events = calculateEventStatistics(obj, regions)
            % required input: the output of a call to findEventRegions
            % calculates statistics about each event and packages it into a
            % cell array of structs
            
            % initialize structure
            events = cell(size(regions,1),1);
            
            % loop through each event
            for i = 1:size(regions,1)
                % grab data
                if diff(regions(i,:))<5e5
                    % just get the whole thing
                    d = obj.sigdata.get(regions(i,:)); % time, current, voltage
                else
                    % have to downsample to save memory
                    d = obj.downsample_pointwise(regions(i,:), 1e5);
                end
                % pad for the std calculation
                pad = round(max(0,min(diff(regions(i,:))/2-5,40)));
                % calculate statistics
                events{i}.current_mean = mean(d(:,2))*obj.in.currentscaling; % current in pA
                events{i}.current_median = median(d(:,2))*obj.in.currentscaling; % current in pA
                events{i}.current_std = std(d(1+pad:end-pad,2))*obj.in.currentscaling; % current in pA
                events{i}.current_range = [min(d(1+pad:end-pad,2)), max(d(1+pad:end-pad,2))]*obj.in.currentscaling; % current in pA
                events{i}.conductance_mean = mean(d(:,2)./d(:,3))*obj.in.currentscaling/obj.in.voltagescaling; % conductance in nS
                events{i}.conductance_median = median(d(:,2)./d(:,3))*obj.in.currentscaling/obj.in.voltagescaling; % conductance in nS
                events{i}.conductance_std = std(d(:,2)./d(:,3))*obj.in.currentscaling/obj.in.voltagescaling; % conductance in nS
                events{i}.voltage = mean(d(:,3))*obj.in.voltagescaling; % voltage in mV
                events{i}.duration = d(end,1)-d(1,1); % duration in seconds
                events{i}.index = regions(i,:);
                events{i}.time = regions(i,:) * obj.sigdata.si;
                % get local open pore value
                open = obj.sigdata.get(max(1,regions(i,1)-100):max(1,regions(i,1)));
                if size(open,1)<2
                    i1 = 0;
                    i2 = 2;
                else
                    open_current = open(:,2)*obj.in.currentscaling; % in pA
                    i2 = find(diff(open_current)>=0,1,'last'); % index of open
                    guessval = open_current(i2);
                    i1 = find(abs(open_current(1:i2))<0.95*abs(guessval),1,'last'); % index of open
                    if isempty(i1)
                        i1 = 0;
                    end
                end
                events{i}.open_pore_current_mean = mean(open_current(i1+1:i2-1)); % pA
                events{i}.open_pore_current_std = std(open_current(i1+1:i2-1)); % pA
                events{i}.open_pore_conductance_mean = mean(open_current(i1+1:i2-1)./open((i1+1):(i2-1),3)); % nS
                events{i}.fractional_block_mean = events{i}.current_mean / events{i}.open_pore_current_mean;
                % check whether the event ended manually (voltage decreased at end)
                d = obj.sigdata.get(regions(i,2) + [1e-4, 1e-3]/obj.sigdata.si); % from 100us after to 1ms after
                v_after = mean(d(:,3))*obj.in.voltagescaling;
                events{i}.voltage_after_event = v_after; % voltage just after event ends
                events{i}.ended_manually = round(v_after) < round(events{i}.voltage); % so did we end it by flipping voltage
                % add file data
                events{i}.file = obj.sigdata.filename;
            end
            
        end
        
        function events = getEvents(obj, varargin)
            % get events struct and return it using other functions
            % optional inputs: trange, voltage, threshold, eventstart, mincond, maxcond
            
            r = obj.findEventRegions(varargin{:});
            events = obj.calculateEventStatistics(r);
            
        end
        
        function logic = getLogic(obj, events, varargin)
            % generate a logical array of size events that says whether
            % each event satisfies the criteria in varargin
            % and the conditions in 'eventlogic' optional input field:
            % eventlogic should be a struct with optional fields that match
            % the names of fields in events, but which are logic operations
            % which must all be true for getLogic to be true
            
            % input handling
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            logic = cellfun(@(x) (isempty(obj.in.files) || any(strcmp(x.file,obj.in.files))) ... % check matching filename
                && (isempty(obj.in.voltage) || any(round(x.voltage/5)*5 == round(obj.in.voltage))), events); % and matching voltage
            
            % get any extra conditionals from 'eventlogic' argument
            if isempty(obj.in.eventlogic)
                return;
            end
            fields = fieldnames(obj.in.eventlogic);
            for i = 1:numel(fields)
                % dynamic field name in eventlogic references function
                % handle which gets passed the argument from events' same
                % field.  this is done for all events.
                condition = cellfun(@(x) obj.in.eventlogic.(fields{i})(x.(fields{i})), events);
                logic = logic & condition; % each time update overall logic
            end
            
        end
        
        function mol = getMolecules(obj, events)
            % get the events which seem to be enzyme-driven molecules
            % this is defined as: not ended manually, I/I_0 < 0.4 && > 0.05
            % and duration > 1 sec
            
            % eventlogic
            eventlogic = struct();
            eventlogic.duration = @(x) x > 1;
            eventlogic.ended_manually = @(x) x == false;
            eventlogic.fractional_block_mean = @(x) x < 0.4 && x > 0.05;
            obj.parsed = false; % let parser know it must re-parse inputs
            
            mol = events(obj.getLogic(events, 'eventlogic', eventlogic));
        end
        
        function events = findEventLevels(obj, events, expectedLevelsPerSecond, falsePositivesPerSecond, varargin)
        % perform level-finding on each event
        % find the levels which are significant based on Kevin Karplus'
        % algorithm. (falsePositivesPerSecond = 1e-4 is typical)
            
            % input handling
            obj.parsed = false;
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            % loop through each event
            for i = 1:numel(events)
                % make sure the right file is loaded
                if ~strcmp(obj.sigdata.filename, events{i}.file)
                    obj.sigdata = SignalData(events{i}.file);
                end
                % grab event data
                pts = round(events{i}.duration * obj.in.filter * 5); % sample at five times filter frequency if possible
                % if too many points, chunk it in the following way: divide
                % a first pass coarse, then divide each of those
                if pts>5e6
                    true_filter = obj.in.filter;
                    obj.in.filter = 100; % hijack this filter setting for now for downsampling
                    data = obj.downsample_pointwise(events{i}.index, min(1e7,5*obj.in.filter*events{i}.duration));
                    data = data(:,1:2);
                    data(:,2) = data(:,2)*obj.in.currentscaling;
                    % level find coarsely based on heavily downsampled data
                    display('Large event... using iterative level finding ========')
                    coarse_levels = karplus_levels(data, 1e-10, 1e-50, 10); % try to find only a few (empirical...)
                    obj.in.filter = true_filter; % return to the real filter setting for fine-grain level finding
                    levels = cell(0);
                    for j = 1:numel(coarse_levels)
                        % level find in each for real, and compile
                        pts = round(coarse_levels{j}.duration / (1/(obj.in.filter*5))); % sample at five times filter frequency if possible
                        data = obj.downsample_pointwise([coarse_levels{j}.start_time coarse_levels{j}.end_time]/obj.sigdata.si, pts);
                        data = data(:,1:2);
                        data(:,2) = data(:,2)*obj.in.currentscaling;
                        newlevs = karplus_levels(data, expectedLevelsPerSecond, falsePositivesPerSecond, obj.in.filter);
                        levels = [levels; newlevs];
                    end
                    display('=====================================================')
                else
                    data = obj.downsample_pointwise(events{i}.index, pts);
                    data = data(:,1:2);
                    data(:,2) = data(:,2)*obj.in.currentscaling;
                    % level find
                    levels = karplus_levels(data, expectedLevelsPerSecond, falsePositivesPerSecond, obj.in.filter);
                end
                % store level data in event struct
                events{i}.levels = levels;
                events{i}.level_finding.expectedLevelsPerSecond = expectedLevelsPerSecond;
                events{i}.level_finding.falsePositivesPerSecond = falsePositivesPerSecond;
                events{i}.level_finding.filter = obj.in.filter;
            end
            
        end
        
        function f = plotEventScatter(obj, events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: title for plot
            % plot a scatter plot of event mean fractional blockages versus
            % durations
            
            % input handling
            obj.parsed = false;
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            % if user used 'eventlogic'
            logic = obj.getLogic(events, 'eventlogic', obj.in.eventlogic);
            
            % if user explicitly entered files and voltages
            logic2 = cellfun(@(x) (isempty(obj.in.files) || any(strcmp(x.file,obj.in.files))) ... % check matching filename
                && (isempty(obj.in.voltage) || any(round(x.voltage/5)*5 == round(obj.in.voltage))), events); % and matching voltage
            
            % limit to these events
            events = events(logic & logic2);
            
            % what to plot
            switch obj.in.eventblockage
                case 'mean'
                    y = cellfun(@(x) x.fractional_block_mean, events);
                case 'first'
                    y = cellfun(@(x) x.levels{find(cellfun(@(y) y.duration>1e-3, x.levels),1,'first')}.current_mean / x.open_pore_current_mean, events);
                case 'last'
                    y = cellfun(@(x) x.levels{end}.current_mean / x.open_pore_current_mean, events);
            end
            
            % plot
            f = figure(obj.in.figure);
            if obj.in.inverted == true
                y = 1-y;
            end
            ended_manually = cellfun(@(x) isfield(x,'ended_manually') && x.ended_manually, events);
            sk23event = cellfun(@(x) x.current_std/x.open_pore_current_mean > 0.05 && x.duration > 1, events); % empirical approximation for real event
            plot(cellfun(@(x) x.duration, events(ended_manually))*1000, y(ended_manually),'x','markersize',5,'color',obj.in.color)
            hold on
            plot(cellfun(@(x) x.duration, events(~ended_manually & sk23event))*1000, y(~ended_manually & sk23event),'.','markersize',25,'color','g')
            plot(cellfun(@(x) x.duration, events(~ended_manually))*1000, y(~ended_manually),'o','markersize',3,'color',obj.in.color)
            set(gca,'xscale','log','fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            ylim([0 1])
            xlim([1e-2 2e5])
            title(obj.in.title)
            xlabel('Duration (ms)')
            if obj.in.inverted == false
                ylabel('I / I_0');
                try
                    annotation('textbox', [0.7 0.9 0 0], 'String', ...
                        char(unique(cellfun(@(x) [x.file(end-27:end-20), '\_', x.file(end-7:end-4)], events, 'uniformoutput', false))), ...
                        'FontSize', 20);
                catch ex
                end
            else
                ylabel('\DeltaI / I_0');
                try
                    annotation('textbox', [0.7 0.25 0 0], 'String', ...
                        char(unique(cellfun(@(x) [x.file(end-27:end-20), '\_', x.file(end-7:end-4)], events, 'uniformoutput', false))), ...
                        'FontSize', 20);
                catch ex
                end
            end
        end
        
        function f = plotInteractiveEventScatter(obj, events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: title for plot, inverted
            % plot a scatter plot of event mean fractional blockages versus
            % durations that you can click on, and will plot individuals
            
            % input handling
            obj.parsed = false;
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            % if user used 'eventlogic'
            logic = obj.getLogic(events, 'eventlogic', obj.in.eventlogic);
            
            % if user explicitly entered files and voltages
            logic2 = cellfun(@(x) (isempty(obj.in.files) || any(strcmp(x.file,obj.in.files))) ... % check matching filename
                && (isempty(obj.in.voltage) || any(round(x.voltage/5)*5 == round(obj.in.voltage))), events); % and matching voltage
            
            % limit to these events
            events = events(logic & logic2);
            
            % what to plot
            switch obj.in.eventblockage
                case 'mean'
                    y = cellfun(@(x) x.fractional_block_mean, events);
                case 'first'
                    y = cellfun(@(x) x.levels{find(cellfun(@(y) y.duration>1e-3, x.levels),1,'first')}.current_mean / x.open_pore_current_mean, events);
                case 'last'
                    y = cellfun(@(x) x.levels{end}.current_mean / x.open_pore_current_mean, events);
            end
            ended_manually = cellfun(@(x) isfield(x,'ended_manually') && x.ended_manually, events);
            duration = cellfun(@(x) x.duration, events)*1000;
            
            % plot
            f = figure(obj.in.figure);
            if obj.in.inverted == true
                y = 1-y;
            end
            
            for i = 1:numel(events)
                if ended_manually(i)
                    dot = plot(duration(i), y(i), 'rx','markersize',5);
                else
                    dot = plot(duration(i), y(i), 'o','markersize',3,'color',obj.in.color);
                end
                set(dot,'ButtonDownFcn',@(~,~) obj.plotEvent(events, i));
                hold on
            end
            set(gca,'xscale','log','fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            ylim([0 1])
            title('Interactive event scatter plot')
            xlabel('Duration (ms)')
            if obj.in.inverted == false
                ylabel('I / I_0');
            else
                ylabel('\DeltaI / I_0');
            end
            box on
            
        end
        
        function f = plotInteractiveEventScatter3(obj, events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: figure, color
            % plot a scatter plot of event mean fractional blockages versus
            % durations that you can click on, and will plot individuals
            % plots I/I0, Irange/I0, and duration
            
            % input handling
            obj.parsed = false;
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            % if user used 'eventlogic'
            logic = obj.getLogic(events, 'eventlogic', obj.in.eventlogic);
            
            % if user explicitly entered files and voltages
            logic2 = cellfun(@(x) (isempty(obj.in.files) || any(strcmp(x.file,obj.in.files))) ... % check matching filename
                && (isempty(obj.in.voltage) || any(round(x.voltage/5)*5 == round(obj.in.voltage))), events); % and matching voltage
            
            % limit to these events
            events = events(logic & logic2);
            
            % plot
            f = figure(obj.in.figure);
            for i = 1:numel(events)
                try
%                     if isfield(events{i},'current_range')
%                         rng = events{i}.current_range;
%                     else
%                         pad = min(diff(events{i}.index)/2-1,20);
%                         d = obj.downsample_pointwise(events{i}.index+[pad,-1*pad],1000);
%                         %stdv = std(d(:,2)*obj.in.currentscaling);
%                         current = d(:,2)*obj.in.currentscaling;
%                         current = current(current < in.threshold * events{i}.open_pore_current_mean);
%                         rng = range(current);
%                     end
                    rng = events{i}.current_std;
                    %rng = diff(events{i}.current_range);
                    if obj.in.inverted == true
                        y = 1 - events{i}.fractional_block_mean;
                    else
                        y = events{i}.fractional_block_mean;
                    end
                    if isfield(events{i},'ended_manually')
                        if events{i}.ended_manually
                            dot = plot3(y, rng/events{i}.open_pore_current_mean, ...
                                events{i}.duration*1000, 'x', 'Color', 'r', 'markersize', 5);
                        else
                            dot = plot3(y, rng/events{i}.open_pore_current_mean, ...
                                events{i}.duration*1000, 'o', 'Color', obj.in.color, 'markersize', 3);
                        end
                    else
                        dot = plot3(y, rng/events{i}.open_pore_current_mean, ...
                            events{i}.duration*1000, 'o', 'Color', obj.in.color, 'markersize', 3);
                    end
%                     if rng/events{i}.open_pore_current_mean > 1
%                         pause();
%                     end
                    set(dot,'ButtonDownFcn',@(~,~) obj.plotEvent(events, i));
                    hold on
                catch ex
                    disp(['Unable to plot event ' num2str(i)])
                end
            end
            set(gca,'zscale','log','fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            xlim([0 1])
            %ylim([0 1])
            title('Interactive event scatter plot')
            zlabel('Duration (ms)')
            if obj.in.inverted == false
                xlabel('I / I_0');
            else
                xlabel('\DeltaI / I_0');
            end
            ylabel('I_s_t_d / I_0')
            grid on
            
        end
        
        function f = plotEvent(obj, events, i, varargin)
            % necessary inputs: events cell struct, output of
            % calculateEventStatistics, and event number of interest i
            % plot event i
            
            % input handling
            if numel(varargin)>1
                obj.parsed = false;
                obj.parseOptionalInputs(varargin{:});
                obj.parsed = true;
            end
            
            pad = max(2e-4/obj.sigdata.si, events{i}.duration/50/obj.sigdata.si); % data points before and after
            
            f = figure;
            if ~strcmp(obj.sigdata.filename,events{i}.file)
                obj.sigdata = SignalData(events{i}.file);
            end
            d = obj.downsample_pointwise(events{i}.index + [-1*pad, pad], 50000); % grab data
            
            % plot either in ms or s depending on scale of event
            timefactor = 1;
            if d(end,1)-d(1,1)<1.5
                timefactor = 1000;
            end
            plot((d(:,1)-events{i}.time(1))*timefactor,d(:,2)*obj.in.currentscaling/events{i}.open_pore_current_mean,'k')
            
            % if there are levels specified, show them
            if isfield(events{i},'levels')
                timing = cell2mat(cellfun(@(x) [x.start_time, x.end_time], events{i}.levels, 'uniformoutput', false));
                means = cellfun(@(x) x.current_mean, events{i}.levels) / events{i}.open_pore_current_mean;
                line((timing'-timing(1,1))*timefactor,(means*[1,1])','LineWidth',2);
            end
            
            set(gca,'fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            title(['Event number ' num2str(i)])
            
            if timefactor==1000
                xlabel('Time (ms)')
            else
                xlabel('Time (s)')
            end
            ylabel('I / I_0')
            ylim([0 1.1])
            xlim([-Inf Inf])
            
            try
                annotation('textbox', [0.7 0.25 0 0], 'String', ...
                    [events{i}.file(end-27:end-20), '\_', events{i}.file(end-7:end-4) ' ' num2str(obj.in.filter) 'Hz'], ...
                    'FontSize', 20);
            catch
                annotation('textbox', [0.8 0.93 0 0], 'String', ...
                    [events{i}.file(1:end-4) ' ' num2str(obj.in.filter) 'Hz'], ... % just take off the suffix
                    'FontSize', 20);
            end
            
        end
        
        function f = plotEventSquiggle(obj, events, i, varargin)
            % f = plotEventSquiggle(events, i, varargin)
            % if an event has levels found, plot the squiggle data: current
            % mean versus level number
            
            % input handling
            if numel(varargin)>1
                obj.parsed = false;
                obj.parseOptionalInputs(varargin{:});
                obj.parsed = true;
            end
            
            f = figure;
            means = cellfun(@(x) x.current_mean, events{i}.levels) / events{i}.open_pore_current_mean;
            plot(1:numel(means),means,'o-','Color',obj.in.color)
            
            title(['Squiggle data: event number ' num2str(i)])
            
            ylabel('I (pA)')
            xlabel('Level number')
            xlim([0 numel(means)+1])
            
            f = obj.polishPlot(f, events{i}, true, true, false);

        end
        
        function events = batch(obj, varargin)
            % batch analysis
            % needed input: files
            % optional input: voltages
            
            % input handling
            obj.parsed = false;
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            display('Batch data analysis:')
            
            % choose time ranges of interest in each file
            for i = 1:numel(obj.in.files)
                if ishandle(1)
                    close(1);
                end
                pv = PoreView(obj.in.files{i});
                drawnow;
                display('Set cursors to region of interest')
                pause();
                timeranges{i} = pv.getCursors();
            end
            
            % do batch analysis
            events = cell(0);
            for i = 1:numel(obj.in.files)
                events = [events; obj.batch_findEventsAndSave(obj.in.files{i}, timeranges{i}, obj.in.voltage, varargin{:})];
            end
        end
        
        function events = inputMetadata(obj, events)
            % inputMetaData(events)
            % allows user to input file metadata which enables easier
            % analysis later based on things like [KCl], temperature,
            % [ATP], [ADPNP], and which pore was used.
            
            inpt = inputdlg({'Pore','Temperature','KCl molarity','ATP molarity','ADPNP molarity','Mg molarity'},'Input metadata',1,{'M2-MspA','25','1.0','0.002','0','0.002'});
            pore = inpt{1};
            values = cellfun(@(x) str2double(x), inpt(2:end));
            for i = 1:numel(events)
                % tack on metadata
                events{i}.pore = pore;
                events{i}.temperature = values(1);
                events{i}.KCl_molarity = values(2);
                events{i}.ATP_molarity = values(3);
                events{i}.ADPNP_molarity = values(4);
                events{i}.Mg_molarity = values(5);
            end
            
        end
        
        function save(obj, events, varargin)
            % save the data
            
            % input handling
            obj.parsed = false;
            obj.parseOptionalInputs(varargin{:});
            obj.parsed = true;
            
            % saving
            if isempty(events)
                display('No events found')
            else
                try
                    if isempty(obj.in.savefile)
                        
                        savefile = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-4) '_events.mat'];
                        % make directory if it doesn't exist
                        if exist(['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-19)],'dir')==0
                            mkdir(['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-19)]);
                        end
                        save(savefile,'events'); % save data
                    else
                        savefile = obj.in.savefile;
                        save(savefile,'events'); % save data
                    end
                    display(['Saved event data in ' savefile])
                catch ex
                    display('Trouble saving to specified directory')
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function parseOptionalInputs(obj, varargin)
            % parse all inputs so all methods can use them easily
            if obj.parsed == true
                return;
            end
            p = inputParser;
            
            % defaults and checks
            defaultFilterFreq = 1000;
            defaultSampleFreq = 5000;
            checkFilterFreq = @(x) all([isnumeric(x), numel(x)==1, x>=10, x<=20000]);
            
            defaultTimeRange = [obj.sigdata.tstart, obj.sigdata.tend]; % start and end of file
            checkTimeRange = @(x) all([all(isnumeric(x)), numel(x)==2, x(1)<x(2), x(1)>=obj.sigdata.tstart, x(2)<=obj.sigdata.tend]);
            
            % getting next figure for default purposes
            f = get(groot,'currentfigure');
            if ~isempty(f)
                defaultFigureNum = f.Number + 1;
            else
                defaultFigureNum = 1;
            end
            
            checkPosNum = @(x) all([all(isnumeric(x)), all(x>=0)]);
            
            checkEventStart = @(x) any([strcmp(x, 'voltagedrop'), strcmp(x, 'currentdrop')]);
            
            % set up the inputs
            addOptional(p, 'filter', defaultFilterFreq, checkFilterFreq); % filter frequency
            addOptional(p, 'sample', defaultSampleFreq, checkFilterFreq); % (down-) sampling frequency
            addOptional(p, 'trange', defaultTimeRange, checkTimeRange); % time range of interest
            addOptional(p, 'mincond', 1.4, checkPosNum); % min open pore conductance
            addOptional(p, 'maxcond', 3, checkPosNum); % max open pore conductance
            addOptional(p, 'minduration', 1e-6, checkPosNum); % min event duration
            addOptional(p, 'threshold', 0.90, checkPosNum); % fraction of open pore event threshold
            addOptional(p, 'voltage', [], @(x) all([all(isnumeric(x)), all(abs(x)>=1)])); % voltage(s) of interest
            addOptional(p, 'eventstart', 'currentdrop', checkEventStart); % what defines start of event
            addOptional(p, 'voltagecheck', @(x) (isnumeric(x) & x>1)); % function to use to check if voltages are okay
            addOptional(p, 'files', [], @(y) all(cellfun(@(x) ischar(x), y))); % cell array of filenames
            addOptional(p, 'title', 'Event scatter plot', @(x) ischar(x)); % title on plots
            addOptional(p, 'figure', defaultFigureNum, @isvalid); % figure to plot things on
            addOptional(p, 'color', 'k', @(x) or(ischar(x),checkPosNum(x))); % color to use in plots
            addOptional(p, 'inverted', false, @(x) islogical(x)); % scatter plots: inverted true plots \Delta I / I_0
            addOptional(p, 'eventlogic', struct(), @(x) isstruct(x)); % logical conditions for selecting events (used by getLogic)
            addOptional(p, 'currentscaling', 1000, checkPosNum); % true current (pA) = recorded current value * currentscaling
            addOptional(p, 'voltagescaling', 1, checkPosNum); % true voltage (mV) = recorded voltage value * voltagescaling
            addOptional(p, 'savefile', [], @(x) ischar(x)); % true voltage (mV) = recorded voltage value * voltagescaling
            addOptional(p, 'eventblockage', 'mean', @(x) any(cellfun(@(y) strcmp(x,y), {'mean','first','last'}))); % what to plot for blockage data
            
            % parse
            parse(p,varargin{:});
            obj.in = p.Results;
        end
        
        function [voltage, current] = getViewData(obj, trange)
            % get the downsampled view data
            % only go back to file if we haven't done this already
            
            if any([numel(obj.voltage_view)==1, numel(obj.current_view)==1, all(trange ~= obj.tr)])
                raw = obj.sigdata.getViewData(trange); % grab downsampled data
                obj.voltage_view = raw(:,3)*obj.in.voltagescaling;
                obj.current_view = raw(:,2)*obj.in.currentscaling;
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
        
        function d = downsample_pointwise(obj, inds, pts)
            %DOWNSAMPLE_POINTWISE does a pointwise downsampling, returning ABOUT 'pts' points
            % downsample data in chunks of 2^20
            
            % if there is no filtered data in sigdata, make one
            if obj.in.filter==10000
                chan = 2;
            else
                if any(strcmp(obj.sigdata.getSignalList(),[num2str(obj.in.filter) 'Hz (IN 0)']))
                    chan = find(strcmp(obj.sigdata.getSignalList(),[num2str(obj.in.filter) 'Hz (IN 0)']),1,'first');
                else
                    f = obj.in.filter; % if you don't do this, there will be some weird reference to this object in sigdata... and it won't work
                    chan = obj.sigdata.addVirtualSignal(@(d) filt_lpb(d,4,f), [num2str(f) 'Hz'], 2); % add filtered channel
                end
            end
            
            d = [];
            numpts = 2^20;
            rep = max(1, round(diff(inds) / pts)); % number of original points per downsampled point
            if (rep < 2)
                d = obj.sigdata.get(inds);
                return;
            end
            chunks = floor(diff(inds)/numpts); % number of full chunks
            if chunks ~= 0
                for i = 1:chunks % do chunks of numpts points
                    fulldata = obj.sigdata.get(inds(1)+(i-1)*numpts:inds(1)+i*numpts-1, [1,chan]); % get chunk
                    d = [d; downsample(fulldata,rep)];
                    clear fulldata
                end
            end
            if inds(2) - (inds(1)+chunks*numpts) >= rep
                fulldata = obj.sigdata.get(inds(1)+chunks*numpts:inds(2), [1,chan]); % the last bit that's not a full chunk
                d = [d; downsample(fulldata,rep)];
            end
        end
        
        function range = getDataRange(obj, inds, channel)
            % load data in chunks and return min and max of entire thing
            numpts = 2^20;
            chunks = floor(diff(inds)/numpts); % number of full chunks
            range = [Inf, -Inf];
            if chunks ~= 0
                for i = 1:chunks % do chunks of numpts points
                    d = obj.sigdata.get(inds(1)+(i-1)*numpts:inds(1)+i*numpts-1); % get chunk
                    range = [min(range(1),min(d(:,channel))), max(range(2),max(d(:,channel)))];
                    clear d
                end
            end
            if diff(inds)-chunks*numpts~=0
                d = obj.sigdata.get(inds(1)+chunks*numpts:inds(2)); % get chunk
                range = [min(range(1),min(d(:,channel))), max(range(2),max(d(:,channel)))];
            end
        end
        
        function events = batch_findEventsAndSave(obj, file, trange, voltage, varargin)
            % to be called from obj.batch with a given file
            disp(file)
            a = analysis(SignalData(file)); % new object
            vs = voltage;
            if isempty(vs)
                if isempty(trange)
                    vs = a.getAppliedVoltages(varargin{:});
                    events = a.getEvents('threshold',0.90,'eventstart','currentdrop','voltage',vs,varargin{:}); % do event finding
                else
                    vs = a.getAppliedVoltages('trange', trange, varargin{:});
                    events = a.getEvents('trange',trange,'threshold',0.90,'eventstart','currentdrop','voltage',vs,varargin{:}); % do event finding
                end
            else
                if isempty(trange)
                    events = a.getEvents('threshold',0.90,'eventstart','currentdrop','voltage',vs,varargin{:}); % do event finding
                else
                    events = a.getEvents('trange',trange,'threshold',0.90,'eventstart','currentdrop','voltage',vs,varargin{:}); % do event finding
                end
            end
            obj.save(events);
        end
        
        function f = polishPlot(obj, f, event, showUpper, showTime, showFilter)
            % do the things all plots get: annotation and sizing
            try
                str = [event.file(end-27:end-20), '\_', event.file(end-7:end-4)];
            catch
                str = event.file(1:end-4);
            end
            if showUpper
                loc = [0.7 0.91 0 0];
            else
                loc = [0.7 0.25 0 0];
            end
            if showTime
                str = [str  ' [' sprintf('%.3g',event.time(1)) '-' sprintf('%.3g',event.time(2)) ']s'];
            end
            if showFilter
                str = [str  ' ' num2str(obj.in.filter) 'Hz'];
            end
            annotation('textbox', loc, 'String', str, 'FontSize', 20);
            set(gca,'fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
        end
        
    end
    
end