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
            in = obj.parseOptionalInputs(varargin{:});
            
            % grab voltage trace
            d = obj.sigdata.getViewData(in.trange);
            voltage = medfilt1(d(:,3),10); % median filtered, rough data
            voltage = voltage(in.voltagecheck(voltage)); % apply the function 'voltagecheck' to see which time ranges have valid voltages
            
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
            in = obj.parseOptionalInputs(varargin{:});
            
            % get the open pore conductance
            [g_m, g_s] = obj.getOpenPoreConductance('mincond', in.mincond, 'maxcond', in.maxcond, 'voltage', in.voltage);
            
            % coarse event finding using a threshold
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
                lowcond = abs(conductance) < abs(g_m) * in.threshold; % regions of conductance below threshold
                dlogic = diff([1; lowcond; 0]);
                startcondition = @(x) abs(x(:,2)*1000./x(:,3)) < abs(g_m) * in.threshold;
                startcondition_inv = @(x) abs(x(:,2)*1000./x(:,3)) > abs(g_m) * in.threshold;
                endcondition = @(x,evtcond) abs(x(:,2)*1000./x(:,3)) < abs(g_m) * in.threshold;
            elseif strcmp(in.eventstart, 'voltagedrop')
                lowvolt = abs(voltage) < abs(V) * in.threshold; % regions of voltage below threshold
                dlogic = diff([1; lowvolt; 0]);
                startcondition = @(x) abs(x(:,3)) < abs(V) * in.threshold;
                startcondition_inv = @(x) abs(x(:,3)) > abs(V) * in.threshold;
                endcondition = @(x,evtcond) abs(x(:,2)*1000./x(:,3)) < abs(g_m) * in.threshold;
            end
            clear voltage;
            [~,possibleStartInds] = findpeaks(double(dlogic > 0),'minpeakheight',0.5,'minpeakdist',5);
            % check to make sure current starts at open pore level
            startsHigh = arrayfun(@(x) nanmax(conductance(x-5:x)) > g_m * in.threshold, possibleStartInds);
            possibleStartInds = possibleStartInds(startsHigh);
            % one end for each start
            conductance = [conductance; nan]; % just so it won't try to go past end
            possibleEndInds = arrayfun(@(x) find(or(conductance(x+2:end) > g_m * in.threshold, isnan(conductance(x+2:end))), 1, 'first'), possibleStartInds) + possibleStartInds + 1;
            % to get rid of the off-by-one errors
            possibleStartInds = possibleStartInds - 1;
            
            % trim out ones that are too short
            too_short = (possibleEndInds-possibleStartInds)*dt < in.minduration;
            possibleStartInds = possibleStartInds(~too_short);
            possibleEndInds = possibleEndInds(~too_short);
            
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
                % find end by next cross of threshold
                if possibleEndInds(i)-possibleStartInds(i) < pad*dt/obj.sigdata.si
                    % end is so close we can search from the start
                    end_inds_thresh = obj.sigdata.findNext(@(x) or(x(:,2)*1000./x(:,3) > g_m * in.threshold, x(:,3)<=0), start_inds(i)+pad);
                    end_inds_thresh = obj.sigdata.findPrev(@(x) or(endcondition(x, evtconductance), x(:,3)<=0), end_inds_thresh+1);
                else
                    % end is far enough that we should work from our
                    % best guess of the end itself
                    end_inds_thresh = obj.sigdata.findNext(@(x) or(x(:,2)*1000./x(:,3) > g_m * in.threshold, x(:,3)<=0), ...
                        (in.trange(1)/dt + max(possibleEndInds(i)-pad, possibleStartInds(i))) * dt/obj.sigdata.si+1);
                    end_inds_thresh = obj.sigdata.findPrev(@(x) or(endcondition(x, evtconductance), x(:,3)<=0), end_inds_thresh+1);
                end
                % end should actually be when current starts to return
                % to open pore
                ending_bit = obj.sigdata.get(max(start_inds(i),end_inds_thresh-20):max(start_inds(i),end_inds_thresh-5),2)*1000;
                end_cond = abs(mean(ending_bit)/V);
                %end_cond_std = std(ending_bit);
                end_inds(i) = obj.sigdata.findPrev(@(x) x(:,2)*1000./x(:,3) < mean([end_cond, end_cond, end_cond, g_m * in.threshold]), end_inds_thresh);
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
            too_short = (end_inds-start_inds)*obj.sigdata.si < in.minduration;
            regions = regions(~too_short,:);
        end
        
        function showEventsInPoreView(obj, pv, r_or_events, how)
            % necessary inputs: PoreView object
            % result of findEventRegions or calculateEventStatistics, how ('all' or 'one')
            % plots the events, one-by-one, in PoreView
            
            if iscell(r_or_events)
                % user passed in events cell struct
                r = cell2mat(cellfun(@(x) x.index, r_or_events, 'uniformoutput', false));
            else
                % user passed in region indices
                r = r_or_events;
            end
            
            pv.clearAxes();
            if strcmp(how,'all')
                for i = 1:size(r,1)
                    % use the y values for whatever signal is currently
                    % plotted
                    y(i,1) = pv.data.get(r(i,1),pv.psigs(1).sigs);
                    y(i,2) = pv.data.get(r(i,2),pv.psigs(1).sigs);
                end
                plot(pv.psigs(1).axes, r(:,1)*pv.data.si,y(:,1),'go')
                plot(pv.psigs(1).axes, r(:,2)*pv.data.si,y(:,2),'rx')
            elseif strcmp(how,'one')
                for i = 1:size(r,1)
                    % use the y values for whatever signal is currently
                    % plotted
                    pv.setView(max(0,sort(r(i,:)+[-100, 100])).*pv.data.si);
                    plot(pv.psigs(1).axes, r(i,1)*pv.data.si, pv.data.get(r(i,1),pv.psigs(1).sigs),'go')
                    plot(pv.psigs(1).axes, r(i,2)*pv.data.si, pv.data.get(r(i,2),pv.psigs(1).sigs),'rx')
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
                pad = round(max(0,min(diff(regions(i,:))/2-5,20)));
                % calculate statistics
                events{i}.current_mean = mean(d(:,2))*1000; % current in pA
                events{i}.current_median = median(d(:,2))*1000; % current in pA
                events{i}.current_std = std(d(1+pad:end-pad,2))*1000; % current in pA
                events{i}.current_range = [min(d(1+pad:end-pad,2)), max(d(1+pad:end-pad,2))]*1000; % current in pA
                events{i}.conductance_mean = mean(d(:,2)./d(:,3))*1000; % conductance in nS
                events{i}.conductance_median = median(d(:,2)./d(:,3))*1000; % conductance in nS
                events{i}.conductance_std = std(d(:,2)./d(:,3))*1000; % conductance in nS
                events{i}.voltage = mean(d(:,3)); % voltage in mV
                events{i}.duration = d(end,1)-d(1,1); % duration in seconds
                events{i}.index = regions(i,:);
                events{i}.time = regions(i,:) * obj.sigdata.si;
                % get local open pore value
                open = obj.sigdata.get((regions(i,1)-100):regions(i,1));
                open_current = open(:,2)*1000; % in pA
                i2 = find(diff(open_current)>=0,1,'last'); % index of open
                guessval = open_current(i2);
                i1 = find(abs(open_current(1:i2))<0.95*abs(guessval),1,'last'); % index of open
                if isempty(i1)
                    i1 = 0;
                end
                events{i}.open_pore_current_mean = mean(open_current(i1+1:i2-1)); % pA
                events{i}.open_pore_current_std = std(open_current(i1+1:i2-1)); % pA
                events{i}.open_pore_conductance_mean = mean(open_current(i1+1:i2-1)./open((i1+1):(i2-1),3)); % nS
                events{i}.fractional_block_mean = events{i}.current_mean / events{i}.open_pore_current_mean;
                % check whether the event ended manually (voltage decreased at end)
                d = obj.sigdata.get(regions(i,2) + [1e-4, 1e-3]/obj.sigdata.si); % from 100us after to 1ms after
                v_after = mean(d(:,3));
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
        
        function f = plotEventScatter(obj, events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: title for plot
            % plot a scatter plot of event mean fractional blockages versus
            % durations
            
            % input handling
            in = obj.parseOptionalInputs(varargin{:});
            
            include_logic = cellfun(@(x) (~isempty(in.files) && any(strcmp(x.file,in.files))) ... % check matching filename
                && (~isempty(in.voltage) && any(round(x.voltage/5)*5 == round(in.voltage))), events); % and matching voltage
            events = events(include_logic); % limit to these events
            
            % plot
            f = figure;
            if in.inverted == true
                y = cellfun(@(x) 1 - x.fractional_block_mean, events);
            else
                y = cellfun(@(x) x.fractional_block_mean, events);
            end
            ended_manually = cellfun(@(x) isfield(x,'ended_manually') && x.ended_manually, events);
            plot(cellfun(@(x) x.duration, events(ended_manually))*1000, y(ended_manually),'ro','markersize',3)
            hold on
            plot(cellfun(@(x) x.duration, events(~ended_manually))*1000, y(~ended_manually),'ko','markersize',3)
            set(gca,'xscale','log','fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            ylim([0 1])
            xlim([1e-2 1e5])
            title(in.title)
            xlabel('Duration (ms)')
            if in.inverted == false
                ylabel('I / I_0');
                annotation('textbox', [0.7 0.9 0 0], 'String', ...
                    char(unique(cellfun(@(x) [x.file(end-27:end-20), '\_', x.file(end-7:end-4)], events, 'uniformoutput', false))), ...
                    'FontSize', 20);
            else
                ylabel('\DeltaI / I_0');
                annotation('textbox', [0.7 0.25 0 0], 'String', ...
                    char(unique(cellfun(@(x) [x.file(end-27:end-20), '\_', x.file(end-7:end-4)], events, 'uniformoutput', false))), ...
                    'FontSize', 20);
            end
        end
        
        function f = plotInteractiveEventScatter(obj, events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: title for plot, inverted
            % plot a scatter plot of event mean fractional blockages versus
            % durations that you can click on, and will plot individuals
            
            % input handling
            in = obj.parseOptionalInputs(varargin{:});
            
            % plot
            f = figure;
            for i = 1:numel(events)
                if (~isempty(in.files) && ~any(strcmp(events{i}.file,in.files))) ... % if user specified filenames and this event doesn't match one
                        || (~isempty(in.voltage) && ~any(round(events{i}.voltage/5)*5 == round(in.voltage))) % or if user specified voltages and this event doesn't match one
                    continue; % skip the rest of this loop iteration
                end
                if in.inverted == true
                    y = 1 - events{i}.fractional_block_mean;
                else
                    y = events{i}.fractional_block_mean;
                end
                if isfield(events{i},'ended_manually') && events{i}.ended_manually
                    dot = plot(events{i}.duration*1000, y,'ro','markersize',3);
                else
                    dot = plot(events{i}.duration*1000, y,'ko','markersize',3);
                end
                set(dot,'ButtonDownFcn',@(~,~) obj.plotEvent(events, i));
                hold on
            end
            set(gca,'xscale','log','fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            ylim([0 1])
            title('Interactive event scatter plot')
            xlabel('Duration (ms)')
            if in.inverted == false
                ylabel('I / I_0');
            else
                ylabel('\DeltaI / I_0');
            end
            
        end
        
        function f = plotInteractiveEventScatter3(obj, events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: figure, color
            % plot a scatter plot of event mean fractional blockages versus
            % durations that you can click on, and will plot individuals
            % plots I/I0, Irange/I0, and duration
            
            % input handling
            in = obj.parseOptionalInputs(varargin{:});
            
            % plot
            f = figure(in.figure);
            for i = 1:numel(events)
                try
                    if isfield(events{i},'current_range')
                        rng = events{i}.current_range;
                    else
                        pad = min(diff(events{i}.index)/2-1,20);
                        d = obj.downsample_pointwise(events{i}.index+[pad,-1*pad],1000);
                        %stdv = std(d(:,2)*1000);
                        current = d(:,2)*1000;
                        current = current(current < in.threshold * events{i}.open_pore_current_mean);
                        rng = range(current);
                    end
                    if in.inverted == true
                        y = 1 - events{i}.fractional_block_mean;
                    else
                        y = events{i}.fractional_block_mean;
                    end
                    dot = plot3(y, rng/events{i}.open_pore_current_mean, ...
                        events{i}.duration*1000, 'o', 'Color', in.color, 'markersize', 3);
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
            if in.inverted == false
                xlabel('I / I_0');
            else
                xlabel('\DeltaI / I_0');
            end
            ylabel('I_s_t_d / I_0')
            grid on
            
        end
        
        function f = plotEvent(obj, events, i)
            % necessary inputs: events cell struct, output of
            % calculateEventStatistics, and event number of interest i
            % plot event i
            
            pad = 100; % data points before and after
            
            f = figure;
            if ~strcmp(obj.sigdata.filename,events{i}.file)
                obj.sigdata = SignalData(events{i}.file);
            end
            d = obj.downsample_pointwise(events{i}.index + [-1*pad, pad], 10000); % grab data
            
            % plot either in ms or s depending on scale of event
            if d(end,1)-d(1,1)<1.5
                plot((d(:,1)-d(pad+1,1))*1000,d(:,2)*1000/events{i}.open_pore_current_mean,'k')
                xlabel('Time (ms)')
            else
                plot((d(:,1)-d(pad+1,1)),d(:,2)*1000/events{i}.open_pore_current_mean,'k')
                xlabel('Time (s)')
            end
            
            set(gca,'fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            title(['Event number ' num2str(i)])
            
            ylabel('I / I_0')
            ylim([0 1.1])
            xlim([-Inf Inf])
            
            annotation('textbox', [0.7 0.25 0 0], 'String', ...
                [events{i}.file(end-27:end-20), '\_', events{i}.file(end-7:end-4)], ...
                'FontSize', 20);
            
        end
        
        function events = batch(obj, varargin)
            % batch analysis
            % needed input: files
            % optional input: voltages
            
            % input handling
            in = obj.parseOptionalInputs(varargin{:});
            
            display('Batch data analysis:')
            
            % choose time ranges of interest in each file
            for i = 1:numel(in.files)
                if ishandle(1)
                    close(1);
                end
                pv = PoreView(in.files{i});
                drawnow;
                display('Set cursors to region of interest')
                pause();
                timeranges{i} = pv.getCursors();
            end
            
            % do batch analysis
            events = cell(0);
            for i = 1:numel(in.files)
                events = [events; obj.batch_findEventsAndSave(in.files{i}, timeranges{i}, in.voltage, varargin{:})];
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
        
        function d = downsample_pointwise(obj, inds, pts)
            %DOWNSAMPLE_POINTWISE does a pointwise downsampling, returning ABOUT 'pts' points
            % downsample data in chunks of 2^20
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
                    fulldata = obj.sigdata.get(inds(1)+(i-1)*numpts:inds(1)+i*numpts-1); % get chunk
                    d = [d; downsample(fulldata,rep)];
                    clear fulldata
                end
            end
            if mod(pts,numpts)~=0
                fulldata = obj.sigdata.get(inds(1)+chunks*numpts:inds(2)); % the last bit that's not a full chunk
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
            
            if isempty(events)
                display('No events found')
            else
                try
                    savefile = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-4) '_events.mat'];
                    save(savefile,'events'); % save data
                    display(['Saved event data in ' savefile])
                catch ex
                    display('Trouble saving to specified directory')
                end
            end
        end
        
    end
    
end