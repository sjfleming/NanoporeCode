classdef analysis < handle
% ANALYSIS
% Class for analysis of nanopore data, or other time-dependent signal data.
% Can be instantiated as an object associated with a particular file.
% Also has static methods which can be called to automate batch file
% analysis, or to plot or further analyze results of previous analyses.
% 
% Stephen Fleming
% 2/23/17
    
    properties
        
        % these properties go with a single file for which an analysis
        % object has been created
        sigdata = []; % SignalData object associated with data file
        in = []; % struct containing parsed inputs
        parsed = false; % has optional user input been parsed?
        metadata = []; % metadata about file
        
    end
    
    methods (Static, Access = public) % can be called from outside an object by specifying files or analysis outputs
        
        % batch file analysis and metadata
        
        function events = batch(varargin)
            % events = batch('files',files_cell_array,'optional...')
            % batch analysis
            % required input: files
            % optional inputs: voltages, event threshold, open pore limits
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            display('Batch data analysis:')
            
            % input metadata
            metadata = analysis.inputMetadata();
            if isempty(metadata)
                display('Operation cancelled.');
                events = cell(0);
                return;
            end
            
            % choose time ranges of interest in each file
            for i = 1:numel(in.files)
                if ishandle(1)
                    close(1);
                end
                pv = PoreView(in.files{i});
                drawnow;
                display('Set cursors to region of interest')
                pause();
                tr = pv.getCursors();
                if isempty(tr) % no selection means whole file
                    tr = [pv.data.tstart, pv.data.tend];
                end
                timeranges{i} = tr;
            end
            
            % do batch analysis
            events = cell(0);
            for i = 1:numel(in.files)
                a = analysis(SignalData(in.files{i})); % create an object instance for this file
                a.in = in; % use parsed inputs
                a.parsed = true; % shortcut unnecessary re-parsing
                a.in.trange = timeranges{i}; % which time range for this file
                a.in.files = in.files{i}; % which file
                a.metadata = metadata;
                events = [events; a.findEventsAndSave()];
            end

        end
        
        function events = attachMetadata(events, varargin)
            % events = attachMetadata(events)
            % allows user to input file metadata which enables easier
            % analysis later based on things like [KCl], temperature,
            % [ATP], [ADPNP], and which pore was used.
            
            if numel(varargin)<1
                md = analysis.inputMetadata();
            else
                md = varargin{1};
            end
            fields = fieldnames(md);
            for i = 1:numel(events)
                for j = 1:numel(fields)
                    % field names in metadata are added to events
                    events{i}.(fields{j}) = md.(fields{j});
                end
            end
            
        end
        
        function files = stephenFileNames(date, nums)
            % return cell array of file names according to SJF Mac computer
            % filesystem of naming
            % files = stephenFileNames('20140210',[31,32,43,44])
            folder = ['/Users/Stephen/Documents/Stephen/Research/Data/Biopore/' date '/'];
            files = cell(0);
            for i = nums
                files{end+1} = [folder date(1:4) '_' date(5:6) '_' date(7:8) '_' sprintf('%04d',i) '.abf'];
            end
        end
        
        function events = loadEvents(files)
            % events = loadEvents(files)
            % loads all the event data in all files passed in cell array
            
            events = cell(0);
            for i = 1:numel(files)
                datafile = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' files{i}(end-27:end-4) '_events.mat'];
                try
                    e = load(datafile);
                catch ex
                    display('Error loading requested event files.');
                    events = cell(0);
                    return;
                end
                events = [events; e.events];
            end
            
        end
        
        % further event analysis: level finding, alignment, etc.
        
        function mol = getMolecules(events)
            % get the events which seem to be enzyme-driven molecules
            % this is defined as: not ended manually, I/I_0 < 0.4 && > 0.05
            % and duration > 1 sec
            
            % eventlogic
            eventlogic = struct();
            eventlogic.duration = @(x) x > 1;
            eventlogic.ended_manually = @(x) x == false;
            eventlogic.fractional_block_mean = @(x) x < 0.4 && x > 0.05;
            
            mol = events(analysis.getLogic(events, 'eventlogic', eventlogic));
        end
        
        function logic = getLogic(events, varargin)
            % generate a logical array of size events that says whether
            % each event satisfies the criteria in varargin
            % and the conditions in 'eventlogic' optional input field:
            % eventlogic should be a struct with optional fields that match
            % the names of fields in events, but which are logic operations
            % which must all be true for getLogic to be true
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            logic = cellfun(@(x) (isempty(in.files) || any(strcmp(x.file,in.files))) ... % check matching filename
                && (isempty(in.voltage) || any(round(x.voltage/5)*5 == round(in.voltage))), events); % and matching voltage
            
            % get any extra conditionals from 'eventlogic' argument
            if isempty(in.eventlogic)
                return;
            end
            fields = fieldnames(in.eventlogic);
            for i = 1:numel(fields)
                % dynamic field name in eventlogic references function
                % handle which gets passed the argument from events' same
                % field.  this is done for all events.
                condition = cellfun(@(x) in.eventlogic.(fields{i})(x.(fields{i})), events);
                logic = logic & condition; % each time update overall logic
            end
            
        end
        
        function totalDuration = getTotalRecordingDuration(events, varargin)
            % totalDuration = getTotalRecordingDuration(events, varargin)
            % returns the total duration of the file recording, at least,
            % from the beginning of the first event to the end of the last
            % event.  sums over all files represented in 'events'.
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            % pare down events according to user input
            logic = analysis.getLogic(events, varargin{:});
            events = events(logic);
            
            % loop through each unique file and add up durations
            totalDuration = 0;
            uniqueFiles = unique(cellfun(@(x) x.file, events, 'UniformOutput', false));
            for i = 1:numel(uniqueFiles)
                elog.file = @(x) strcmp(x, uniqueFiles{i});
                logic = analysis.getLogic(events, 'eventlogic', elog);
                t1 = min(cellfun(@(x) x.time(1), events(logic)));
                t2 = max(cellfun(@(x) x.time(2), events(logic)));
                totalDuration = totalDuration + (t2-t1);
            end
            
        end
        
        function events = findEventLevels(events, expectedLevelsPerSecond, falsePositivesPerSecond, varargin)
        % perform level-finding on each event
        % find the levels which are significant based on Kevin Karplus'
        % algorithm. (falsePositivesPerSecond = 1e-4 is typical)
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            % create an object
            a = analysis(SignalData(events{1}.file));
            a.in = in;
            a.parsed = true;
            
            % loop through each event
            for i = 1:numel(events)
                % make sure the right file is loaded
                if ~strcmp(a.sigdata.filename, events{i}.file)
                    a = analysis(SignalData(events{i}.file)); % new object
                    a.in = in;
                    a.parsed = true;
                end
                % grab event data
                pts = round(events{i}.duration * in.filter * 5); % sample at five times filter frequency if possible
                % if too many points, chunk it in the following way: divide
                % a first pass coarse, then divide each of those
                if pts>5e6
                    true_filter = in.filter;
                    a.in.filter = 100; % hijack this filter setting for now for downsampling
                    data = a.downsample_pointwise(events{i}.index, min(1e7,5*a.in.filter*events{i}.duration));
                    data = data(:,1:2);
                    data(:,2) = data(:,2)*in.currentscaling;
                    % level find coarsely based on heavily downsampled data
                    display('Large event... using iterative level finding ========')
                    coarse_levels = karplus_levels(data, 1e-10, 1e-50, 10); % try to find only a few (empirical...)
                    a.in.filter = true_filter; % return to the real filter setting for fine-grain level finding
                    levels = cell(0);
                    for j = 1:numel(coarse_levels)
                        % level find in each for real, and compile
                        pts = round(coarse_levels{j}.duration / (1/(in.filter*5))); % sample at five times filter frequency if possible
                        data = obj.downsample_pointwise([coarse_levels{j}.start_time coarse_levels{j}.end_time]/obj.sigdata.si, pts);
                        data = data(:,1:2);
                        data(:,2) = data(:,2)*in.currentscaling;
                        newlevs = karplus_levels(data, expectedLevelsPerSecond, falsePositivesPerSecond, in.filter);
                        levels = [levels; newlevs];
                    end
                    display('=====================================================')
                else
                    data = a.downsample_pointwise(events{i}.index, pts);
                    data = data(:,1:2);
                    data(:,2) = data(:,2)*in.currentscaling;
                    % level find
                    levels = karplus_levels(data, expectedLevelsPerSecond, falsePositivesPerSecond, in.filter);
                end
                % store level data in event struct
                events{i}.levels = levels;
                events{i}.level_finding.expectedLevelsPerSecond = expectedLevelsPerSecond;
                events{i}.level_finding.falsePositivesPerSecond = falsePositivesPerSecond;
                events{i}.level_finding.filter = in.filter;
            end
            
        end
        
        function events = doIterativeScalingAlignment(events)
            % iterates between model scaling and level alignment until
            % convergence is reached.  helps improve model scaling.
            
            % sequence
            seq = 'RRRRRTTTTTRRRRGGTTGTTTCTGTTGGTGCTGATATTGCGGCGTCTGCTTGGGTGTTTAACCT'; % SK23
            
            for i = 1:numel(events)
                
                % get the initial predicted levels from oxford
                events{i}.sequence = seq;
                levs = cellfun(@(x) x.current_mean, events{i}.levels);
                hicut = 0.6 * abs(events{i}.open_pore_current_mean);
                lowcut = 0.05 * abs(events{i}.open_pore_current_mean);
                [model_levels, model_levels_std] = ...
                    get_model_levels_oxford(events{i}.sequence, levs(levs>lowcut & levs<hicut), ...
                    abs(events{i}.open_pore_current_mean), abs(events{i}.voltage), events{i}.temperature);
%                 [model_levels, model_levels_std,~,~] = ...
%                     get_model_levels_M2(events{i}.sequence, levs(levs>lowcut & levs<hicut));
                
                % save the initial scaling
                events{i}.model_levels = cell(numel(model_levels),1);
                for j = 1:numel(model_levels)
                    events{i}.model_levels{j}.mean = model_levels(j);
                    events{i}.model_levels{j}.stdev = model_levels_std(j);
                end
                
                % iterate alignment and scaling until convergence is achieved
                stop_criterion = 0.05;
                lsqdist = [];
                dfit = 1;
                while dfit > stop_criterion && numel(lsqdist) < 20
                    % do a level alignment
                    events{i} = obj.doLevelAlignment(events{i});
                    % calculate least-squares distance per measured level
                    logic = ~isnan(events{i}.level_alignment.model_levels_measured_mean_currents); % these levels are not missing
                    lsqdist(end+1) = sum((events{i}.level_alignment.model_levels_measured_mean_currents(logic) ...
                        - cellfun(@(x) x.mean, events{i}.model_levels(logic))).^2) / sum(logic);
                    if numel(lsqdist)>1
                        dfit = abs(lsqdist(end)-lsqdist(end-1));
                    end
                    if sum(logic)<2
                        display('Cannot do iterative fitting, too far off.')
                        break;
                    end
                    % do a least-squares fit
                    f = fit(events{i}.level_alignment.model_levels_measured_mean_currents(logic), ...
                        cellfun(@(x) x.mean, events{i}.model_levels(logic)), ...
                        'poly1','Weights', ...
                        events{i}.level_alignment.model_levels_measured_total_duration(logic), ...
                        'Robust','bisquare');
                    % invert to get the scale and offset corrections and update the
                    % predicted levels
                    for j = 1:numel(events{i}.model_levels)
                        events{i}.model_levels{j}.mean = (events{i}.model_levels{j}.mean-f.p2)/f.p1;
                        events{i}.model_levels{j}.stdev = events{i}.model_levels{j}.stdev/f.p1;
                    end
                end
                
                if numel(lsqdist)==20
                    display('Iterative scaling alignment warning: convergence not reached in 20 iterations.')
                end
                
            end
            
        end
        
        function event = doLevelAlignment(event)
            % align the levels to the predicted model levels
            % error if we are missing something
            if ~isfield(event,'levels')
                display('Error: could not align levels.  No levels extracted from the data!');
                return;
            end
            if ~isfield(event,'model_levels')
                display('Error: could not align levels.  No predicted sequence levels!');
                return;
            end
            event.level_alignment = struct(); % clear any previous alignment
            % get reasonable levels
            [mod_inds, mod_type, lvl_accum, P, ks] = align_fb(cellfun(@(x) x.mean, event.model_levels), ...
                cellfun(@(x) x.stdev, event.model_levels), cellfun(@(x) x.current_mean, event.levels), ...
                cellfun(@(x) x.duration, event.levels), 0.18*abs(event.open_pore_current_mean));
            event.level_alignment.model_level_assignment = mod_inds;
            event.level_alignment.level_type = mod_type; % 1 is normal, 2 is noise, 3 is deep block
            event.level_alignment.model_levels_measured_mean_currents = lvl_accum; % mean of level currents for each level assigned to a given model level
            event.level_alignment.model_levels_measured_total_duration = accumarray(mod_inds(mod_type~=2), ...
                cellfun(@(x) x.duration, event.levels(mod_type~=2)), ...
                size(event.model_levels),@sum,nan); % total time in each level
            event.level_alignment.P = P; % probabilities in the alignment matrix
            event.level_alignment.ks = ks; % state index in the alignment matrix
            
            % also store the level means and timing after alignment, by
            % which i mean: combine adjacent "stay" levels
            lmean = [];
            lmed = [];
            lstd = [];
            ltime = [];
            tmp = cellfun(@(x) x.current_mean, event.levels);
            level_means_no_noise = tmp(mod_type==1);
            tmp = cellfun(@(x) x.current_median, event.levels);
            level_medians_no_noise = tmp(mod_type==1);
            tmp = cellfun(@(x) x.current_std, event.levels);
            level_stds_no_noise = tmp(mod_type==1);
            tmp = cell2mat(cellfun(@(x) [x.start_time x.end_time], event.levels, 'uniformoutput', false));
            level_timing_no_noise = tmp(mod_type==1,:);
            mod_inds_no_noise = mod_inds(mod_type==1);
            % go through each level
            i = 1;
            while i <= numel(mod_inds_no_noise)
                % find the region of contiguous stay levels
                last_contiguous_same_ind = find(mod_inds_no_noise(i+1:end) ~= mod_inds_no_noise(i),1,'first') + i - 1;
                if isempty(last_contiguous_same_ind)
                    last_contiguous_same_ind = i;
                end
                region = (1:numel(mod_inds_no_noise)) >= i & (1:numel(mod_inds_no_noise)) <= last_contiguous_same_ind;
                % compute stats
                durations = diff(level_timing_no_noise(region,:),1,2);
                lmean(end+1) = sum(durations .* level_means_no_noise(region)) / sum(durations);
                lmed(end+1) = sum(durations .* level_medians_no_noise(region)) / sum(durations);
                lstd(end+1) = sum(durations .* level_stds_no_noise(region)) / sum(durations);
                if size(ltime,1)<2
                    ltime(end+1,:) = [level_timing_no_noise(i,1), level_timing_no_noise(last_contiguous_same_ind,2)];
                else
                    ltime(end+1,:) = [ltime(end,2), level_timing_no_noise(last_contiguous_same_ind,2)];
                end
                % update index location
                i = last_contiguous_same_ind + 1;
            end
            % save the values in obj.alignment
            event.level_alignment.level_means = lmean';
            event.level_alignment.level_medians = lmed';
            event.level_alignment.level_stds = lstd';
            event.level_alignment.level_timing = ltime - ltime(1,1); % zero at beginning of molecule
            
        end
        
        % plotting of all kinds
        
        function showEventsInPoreView(pv, events, how)
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
                        means = cellfun(@(x) x.current_mean, events{i}.levels) / 1000; % pA back to initial data value
                        plot(pv.psigs(1).axes, timing', (means*[1,1])', '-', 'LineWidth', 4, 'Color', 'k');
                    end
                    pause();
                end
            end
        end
        
        function f = plotEventScatter(events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: title for plot
            % plot a scatter plot of event mean fractional blockages versus
            % durations
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            % get logical array for which events are specified
            logic = analysis.getLogic(events, varargin{:});
            
            % limit to these events
            events = events(logic);
            if isempty(events)
                display('No events to display!')
                return;
            end
            
            % what to plot
            switch in.eventblockage
                case 'mean'
                    y = cellfun(@(x) x.fractional_block_mean, events);
                case 'first'
                    y = cellfun(@(x) x.levels{find(cellfun(@(y) y.duration>1e-3, x.levels),1,'first')}.current_mean / x.open_pore_current_mean, events);
                case 'last'
                    y = cellfun(@(x) x.levels{end}.current_mean / x.open_pore_current_mean, events);
            end
            
            % plot
            f = figure(in.figure);
            if in.inverted == true
                y = 1-y;
            end
            ended_manually = cellfun(@(x) isfield(x,'ended_manually') && x.ended_manually, events);
            %sk23event = cellfun(@(x) x.current_std/x.open_pore_current_mean > 0.05 && x.duration > 1, events); % empirical approximation for real event
            plot(cellfun(@(x) x.duration, events(ended_manually))*1000, y(ended_manually),'o','markersize',6,'color',in.color)
            hold on
            %plot(cellfun(@(x) x.duration, events(~ended_manually & sk23event))*1000, y(~ended_manually & sk23event),'.','markersize',25,'color','g')
            plot(cellfun(@(x) x.duration, events(~ended_manually))*1000, y(~ended_manually),'o','markersize',1.5,'color',in.color)
            set(gca,'xscale','log','fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            ylim([0 1])
            xlim([1e-2 2e5])
            % title(in.title)
            title([in.title ', ' num2str(sum(logic)) ' events in ' num2str(round(analysis.getTotalRecordingDuration(events, varargin{:}))) ' sec'])
            xlabel('Duration (ms)')
            if in.inverted == false
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
        
        function f = plotInteractiveEventScatter(events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: title for plot, inverted
            % plot a scatter plot of event mean fractional blockages versus
            % durations that you can click on, and will plot individuals
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            % if user used 'eventlogic'
            logic = analysis.getLogic(events, varargin{:});
            
            % limit to these events
            events = events(logic);
            if isempty(events)
                display('No events to display!')
                return;
            end
            
            % what to plot
            switch in.eventblockage
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
            f = figure(in.figure);
            if in.inverted == true
                y = 1-y;
            end
            
            sd = SignalData(events{1}.file);
            a = analysis(sd);
            for i = 1:numel(events)
                % load SignalData
                if ~strcmp(sd.filename,events{i}.file)
                    sd = SignalData(events{i}.file);
                    a = analysis(sd);
                end
                % plot
                if ended_manually(i)
                    dot = plot(duration(i), y(i), 'rx','markersize',5);
                else
                    dot = plot(duration(i), y(i), 'o','markersize',3,'color',in.color);
                end
                set(dot,'ButtonDownFcn',@(~,~) a.plotSingleEvent(events{i},varargin{:}));
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
            box on
            
            function keyFn(e, events, originalfig)
            % keyboard callback
                % do this even if we don't have data loaded
                if strcmp(e.Character,'g')
                    % allow user to select a rectangular area            
                    try
                        rect = getrect(gcf);
                        minduration = (rect(1))*1e-3;
                        maxduration = (rect(1)+rect(3))*1e-3;
                        minblock = rect(2);
                        maxblock = rect(2)+rect(4);
                        evtlogic.duration = @(x) x>minduration && x<maxduration;
                        evtlogic.fractional_block_mean = @(x) x>minblock && x<maxblock;
                        evtsSelected = events(analysis.getLogic(events, 'eventlogic', evtlogic));
                        cmap = get(gca,'colororder');
                        c = cmap(randi(size(cmap,1)),:);
                        % change the color of the events selected
                        for num = 1:numel(originalfig.Children.Children)
                            if any(cellfun(@(x) strcmp('Marker',x), fieldnames(originalfig.Children.Children(num)))) % is it a data point obj
                                mx = originalfig.Children.Children(num).XData;
                                my = originalfig.Children.Children(num).YData;
                                if evtlogic.duration(mx*1e-3) && evtlogic.fractional_block_mean(my)
                                    set(originalfig.Children.Children(num),'Color',c);
                                end
                            end
                        end
                        % draw box around them
                        rectangle('Position',rect,'EdgeColor',c);
                        % do the scatter plot overlay of selected events
                        analysis.plotEventsOverlaid(evtsSelected);
                    catch
                        disp('Problem with rectangle selection.  Try again.')
                    end
                    
                end
            end
            set(f,'WindowKeyPressFcn',@(~,e) keyFn(e, events, f));
            
        end
        
        function f = plotInteractiveEventScatter3(events, varargin)
            % necessary input: events cell struct, output of
            % calculateEventStatistics
            % optional input: figure, color
            % plot a scatter plot of event mean fractional blockages versus
            % durations that you can click on, and will plot individuals
            % plots I/I0, Irange/I0, and duration
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            % if user used 'eventlogic'
            logic = analysis.getLogic(events, varargin{:});
            
            % limit to these events
            events = events(logic);
            
            % plot
            f = figure(in.figure);
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
                                events{i}.duration*1000, 'o', 'Color', in.color, 'markersize', 3);
                        end
                    else
                        dot = plot3(y, rng/events{i}.open_pore_current_mean, ...
                            events{i}.duration*1000, 'o', 'Color', in.color, 'markersize', 3);
                    end
%                     if rng/events{i}.open_pore_current_mean > 1
%                         pause();
%                     end
                    set(dot,'ButtonDownFcn',@(~,~) obj.plotSingleEvent(events{i},varargin{:}));
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
        
        function f = plotEvent(event, varargin)
            % f = plotEvent(event, 'optional_inputs'...)
            % necessary inputs: event cell struct
            % optional: filter, color, current scaling, sampling, etc.
            % returns figure handle
            
            % create an analysis object for obtaining the event data
            a = analysis(SignalData(event.file));
            f = a.plotSingleEvent(event, varargin{:});
            
        end
        
        function f = plotEventsOverlaid(events, varargin)
            % f = plotEventsOverlaid(events, 'optional_inputs'...)
            % necessary inputs: event cell struct
            % optional: filter, color, current scaling, sampling, etc.
            % returns figure handle
            
            % create an analysis object for obtaining the event data
            f = figure();
            a = analysis(SignalData(events{1}.file));
            for i = 1:numel(events)
                if ~strcmp(a.sigdata.filename,events{i}.file)
                    a = analysis(SignalData(events{i}.file));
                end
                f = a.plotSingleEventScatter(events{i}, 'figure', f, varargin{:});
            end
            
        end
        
        function f = plotEventSquiggle(event, varargin)
            % f = plotEventSquiggle(event, varargin)
            % if an event has levels found, plot the squiggle data: current
            % mean versus level number.  return figure handle f.
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            f = in.figure;
            means = cellfun(@(x) x.current_mean, event.levels) / event.open_pore_current_mean;
            plot(1:numel(means),means,'o-','Color',in.color)
            
            if isempty(in.title)
                title('Squiggle data')
            else
                title(in.title)
            end
            
            ylabel('I (pA)')
            xlabel('Level number')
            xlim([0 numel(means)+1])
            
            f = analysis.finishPlot(f, event, true, true, false);

        end
        
        function f = plotAlignment(event, varargin)
            % plot the alignment information
            
            % check for alignment
            if ~isfield(event.level_alignment,'model_level_assignment')
                %obj.do_iterative_scaling_alignment;
                display('Alignment not completed!')
                return;
            end
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            f = in.figure();
            
            mod_inds = event.level_alignment.model_level_assignment;
            mod_type = event.level_alignment.level_type; % 1 is normal, 2 is noise, 3 is deep block
            lvl_accum = event.level_alignment.model_levels_measured_mean_currents; % mean of level currents for each level assigned to a given model level
            lvls = cellfun(@(x) x.current_mean, event.levels);
            
            subplot(3,1,1);
            plot(cellfun(@(x) x.mean, event.model_levels),'o-','LineWidth',2)
            hold on
            plot(lvl_accum,'o-','LineWidth',2);
            legend('Model','Data')
            ylabel('Current (pA)')
            xlabel('Model level')
            xlim([0 find(~isnan(lvl_accum),1,'last')+1])
            set(gca,'FontSize',20)
            title('Best fit of data to model')
            
            subplot(3,1,2);
            plot(abs(lvls))
            for i=1:numel(lvls)
                text(i,abs(lvls(i)),num2str(mod_inds(i)),'FontSize',14);
            end
            hold on
            xx = 1:numel(event.levels);
            plot(xx(mod_type==2),cellfun(@(x) x.current_mean, event.levels(mod_type==2)),'rx','MarkerSize',10)
            plot(xx(mod_type==3),cellfun(@(x) x.current_mean, event.levels(mod_type==3)),'go','MarkerSize',10)
            ylabel('Current (pA)')
            xlabel('Measured level')
            xlim([0 numel(lvls)+1])
            set(gca,'FontSize',20)
            title('Matching each measured level to a model state')
            
            subplot(3,1,3);
            imagesc(1.03.^(event.level_alignment.P/2)) % scales so image shows up well
            hold on
            plot(1:numel(lvls),event.level_alignment.ks,'r','LineWidth',3);
            title('Probability Matrix')
            ylabel('State')
            xlabel('Level Number')
            set(gca,'FontSize',20)
            
        end
        
        % saving data
        
        function save(events, varargin)
            % save the data
            
            % input handling
            in = analysis.parseInputs(varargin{:});
            
            % saving
            if isempty(events)
                display('No events found')
            else
                try
                    if isempty(in.savefile)
                        savefile = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-4) '_events.mat'];
                        % make directory if it doesn't exist
                        if exist(['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-19)],'dir')==0
                            mkdir(['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' events{1}.file(end-27:end-19)]);
                        end
                    else
                        savefile = in.savefile;
                    end
                    save(savefile,'events'); % save data
                    display(['Saved event data in ' savefile])
                catch ex
                    display(['Trouble saving to specified directory ' savefile])
                end
            end
        end
        
    end
    
    methods (Static, Access = private) % called outside an object, but not by a user directly
        
        function in = parseInputs(varargin)
            % parse all inputs so all methods can use them easily
            p = inputParser;
            
            % defaults and checks
            defaultFilterFreq = 10000;
            defaultSampleFreq = 50000;
            checkFilterFreq = @(x) all([isnumeric(x), numel(x)==1, x>=10, x<=20000]);
            
            defaultTimeRange = [];
            checkTimeRange = @(x) all([all(isnumeric(x)), numel(x)==2, x(1)<x(2)]);
            
            checkPosNum = @(x) all([all(isnumeric(x)), all(x>=0)]);
            
            checkEventStart = @(x) any([strcmp(x, 'voltagedrop'), strcmp(x, 'currentdrop')]);
            
            % getting next figure for default purposes
            a = 1;
            f = get(groot,'currentfigure');
            while ~isempty(f) && ishandle(a)
                a = a+1;
            end
            defaultFigureNum = a;
            
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
            addOptional(p, 'eventlogic', struct(), @(x) isstruct(x)); % logical conditions for selecting events (used by getLogic)
            addOptional(p, 'currentscaling', 1000, checkPosNum); % true current (pA) = recorded current value * currentscaling
            addOptional(p, 'voltagescaling', 1, checkPosNum); % true voltage (mV) = recorded voltage value * voltagescaling
            addOptional(p, 'savefile', '', @(x) ischar(x)); % true voltage (mV) = recorded voltage value * voltagescaling
            addOptional(p, 'title', 'Scatter plot', @(x) ischar(x)); % title on plots
            addOptional(p, 'figure', defaultFigureNum, @isvalid); % figure to plot things on
            addOptional(p, 'color', 'k', @(x) or(ischar(x),checkPosNum(x))); % color to use in plots
            addOptional(p, 'inverted', false, @(x) islogical(x)); % scatter plots: inverted true plots \Delta I / I_0
            addOptional(p, 'eventblockage', 'mean', @(x) any(cellfun(@(y) strcmp(x,y), {'mean','first','last'}))); % what to plot for blockage data
            
            % parse
            parse(p,varargin{:});
            in = p.Results;
        end
        
        function metadata = inputMetadata()
            % inputMetaData()
            % prompts user to input file metadata which enables easier
            % analysis later based on things like [KCl], temperature,
            % [ATP], [ADPNP], and which pore was used.
            % metadata is added to events that are found later.
            
            inpt = inputdlg({'Pore','Temperature','KCl molarity','ATP molarity','ADPNP molarity','Mg molarity','Analyte'},'Input metadata',1,{'M2-MspA','25','1.0','0.002','0','0.002','E5/SK23-24'});
            values = inpt;
            
            for i = 1:numel(inpt)
                if ~isnan(str2double(inpt{i}))
                    values{i} = str2double(inpt{i});
                end
            end
            
            if isempty(values)
                % user canceled operation
                metadata = [];
                return;
            end
            
            metadata.pore = values{1};
            metadata.temperature = values{2};
            metadata.KCl_molarity = values{3};
            metadata.ATP_molarity = values{4};
            metadata.ADPNP_molarity = values{5};
            metadata.Mg_molarity = values{6};
            metadata.analyte = values{7};
            
        end
        
        function f = finishPlot(f, event, showUpper, showTime, showFilter, filter)
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
                str = [str  ' ' num2str(filter) 'Hz'];
            end
            annotation('textbox', loc, 'String', str, 'FontSize', 20);
            set(gca,'fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
        end
        
    end
    
    methods % can be called directly by the user on a particular analysis object
        
        % constructor
        function obj = analysis(sigdata)
            % initialize the analysis object by linking it to one 
            % SignalData object
            
            obj.sigdata = sigdata;
            
        end
        
        function events = getEvents(obj, varargin)
            % get events struct and return it using other functions
            % optional inputs: trange, voltage, threshold, eventstart, mincond, maxcond
            
            % input handling
            obj.parseObjectInputs(varargin{:});
            
            r = obj.findEventRegions();
            events = obj.calculateEventStatistics(r);
            
        end
        
        function [cond, cond_std] = getOpenPoreConductance(obj, varargin)
            % return open pore conductance
            % based on a histogram
            % optional inputs: trange, voltage, mincond, maxcond
            
            % input handling
            obj.parseObjectInputs(varargin{:});
            
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
            obj.parseObjectInputs(varargin{:});
            
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
        
        function f = plotSingleEvent(obj, event, varargin)
            % f = plotSingleEvent(event, 'optional_inputs'...)
            % necessary inputs: event cell struct
            % optional: filter, color, current scaling, sampling, etc.
            % returns figure handle
            
            % input handling
            obj.parseObjectInputs(varargin{:});
            
            f = figure(obj.in.figure);
            
            pad = max(2e-4/obj.sigdata.si, event.duration/50/obj.sigdata.si); % data points before and after
            
            d = obj.downsample_pointwise(event.index + [-1*pad, pad], 50000); % grab data
            
            % plot either in ms or s depending on scale of event
            timefactor = 1;
            if d(end,1)-d(1,1)<1.5
                timefactor = 1000;
            end
            plot((d(:,1)-event.time(1))*timefactor,d(:,2)*obj.in.currentscaling/event.open_pore_current_mean,'k')
            
            % if there are levels specified, show them
            if isfield(event,'levels')
                timing = cell2mat(cellfun(@(x) [x.start_time, x.end_time], event.levels, 'uniformoutput', false));
                means = cellfun(@(x) x.current_mean, event.levels) / event.open_pore_current_mean;
                line((timing'-timing(1,1))*timefactor,(means*[1,1])','LineWidth',2);
            end
            
            title('Event trace')
            
            if timefactor==1000
                xlabel('Time (ms)')
            else
                xlabel('Time (s)')
            end
            ylabel('I / I_0')
            ylim([0 1.1])
            xlim([-Inf Inf])
            
            analysis.finishPlot(f, event, true, true, true, obj.in.filter);
            
        end
        
        function f = plotSingleEventScatter(obj, event, varargin)
            % f = plotSingleEvent(event, 'optional_inputs'...)
            % necessary inputs: event cell struct
            % optional: filter, color, current scaling, sampling, etc.
            % returns figure handle
            
            % input handling
            obj.parseObjectInputs(varargin{:});
            
            f = figure(obj.in.figure);
            
            pad = max(2e-4/obj.sigdata.si, event.duration/50/obj.sigdata.si); % data points before and after
            
            d = obj.downsample_pointwise(event.index + [-1*pad, pad], 50000); % grab data
            
            % plot either in ms or s depending on scale of event
            timefactor = 1;
            if d(end,1)-d(1,1)<1.5
                timefactor = 1000;
            end
            plot((d(:,1)-event.time(2))*timefactor,d(:,2)*obj.in.currentscaling/event.open_pore_current_mean,'k.')
            
            % if there are levels specified, show them
            if isfield(event,'levels')
                timing = cell2mat(cellfun(@(x) [x.start_time, x.end_time], event.levels, 'uniformoutput', false));
                means = cellfun(@(x) x.current_mean, event.levels) / event.open_pore_current_mean;
                line((timing'-timing(1,1))*timefactor,(means*[1,1])','LineWidth',2);
            end
            
            title('Event trace')
            
            if timefactor==1000
                xlabel('Time (ms)')
            else
                xlabel('Time (s)')
            end
            ylabel('I / I_0')
            ylim([0 1.1])
            xlim([-Inf Inf])
            hold on
            set(gca,'fontsize',18,'LooseInset',[0 0 0 0],'OuterPosition',[0 0 0.99 1])
            
        end
        
    end
    
    methods (Access = private) % only called by class methods
        
        function parseObjectInputs(obj, varargin)
            % parse all inputs so all methods can use them easily
            if obj.parsed == true
                % getting next figure for default purposes
                if ~(any(cellfun(@(x) strcmp('figure',x), varargin)))
                    % figure not specified
                    % so increment to next figure
                    a = 1;
                    f = get(groot,'currentfigure');
                    while ~isempty(f) && ishandle(a)
                        a = a+1;
                    end
                    obj.in.figure = a;
                end
                return;
            end
            obj.in = analysis.parseInputs(varargin{:});
            obj.parsed = true;
        end
        
        function [voltage, current] = getViewData(obj, trange)
            % get the downsampled view data
            
            raw = obj.sigdata.getViewData(trange); % grab downsampled data
            voltage = raw(:,3)*obj.in.voltagescaling;
            current = raw(:,2)*obj.in.currentscaling;
            
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
        
        function regions = findEventRegions(obj)
            % find events in data specified by user parameters
            % optional inputs: trange, voltage, threshold, eventstart, mincond, maxcond
            
            % get the open pore conductance
            [g_m, g_s] = obj.getOpenPoreConductance();
            
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
                [voltage_view,~] = obj.getViewData(obj.in.trange);
                V = mode(round(voltage_view(max(1,possibleStartInds(i)):min(numel(voltage_view),possibleEndInds(i))))); % estimate of this event's voltage
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
                d = obj.sigdata.get(inds, [1,chan,3]);
                return;
            end
            chunks = floor(diff(inds)/numpts); % number of full chunks
            if chunks ~= 0
                for i = 1:chunks % do chunks of numpts points
                    fulldata = obj.sigdata.get(inds(1)+(i-1)*numpts:inds(1)+i*numpts-1, [1,chan,3]); % get chunk
                    d = [d; downsample(fulldata,rep)];
                    clear fulldata
                end
            end
            if inds(2) - (inds(1)+chunks*numpts) >= rep
                fulldata = obj.sigdata.get(inds(1)+chunks*numpts:inds(2), [1,chan,3]); % the last bit that's not a full chunk
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
        
        function events = findEventsAndSave(obj)
            % called by static method 'batch' on a particular analysis
            % object
            disp(obj.in.files)
            if isempty(obj.in.voltage)
                obj.in.voltage = a.getAppliedVoltages();
            end
            events = obj.getEvents(obj.in.trange); % do event finding
            events = obj.attachMetadata(events, obj.metadata); % metadata
            obj.save(events,'savefile',obj.in.savefile);
        end
        
    end
    
end