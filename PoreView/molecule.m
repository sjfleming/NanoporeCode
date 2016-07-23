classdef molecule < handle & matlab.mixin.SetGet
    %MOLECULE object to aid in analysis of nanopore sequencing data
    % Stephen Fleming
    
    properties (SetAccess = public, GetAccess = public)
        start_ind = nan;
        end_ind = nan;
        start_time = nan;
        end_time = nan;
        ended_manually = false;
        start_file = [];
        end_file = [];
        sequence = [];
        predicted_levels = nan;
        predicted_levels_stdev = nan;
        level_means = nan;
        level_medians = nan;
        level_stds = nan;
        level_timing = nan;
        level_finding_params = struct();
        level_alignment = struct();
        pulses = nan;
        voltage = nan;
        temp = nan;
        KCl_molarity = nan;
        ADPNP_molarity = nan;
        open_pore_current = nan;
    end
    
    methods
        
        % constructor, gets passed in a struct
        function obj = molecule(event)
            
            % check if we were passed in a file
            if ~isstruct(event)
                load(event);
                if ~isstruct(molecule_data)
                    display('Error: Attempted to load a non-molecule data file.')
                    return;
                end
                event = molecule_data;
            end
            
            % load in the data
            if isfield(event,'start_ind')
                obj.start_ind = event.start_ind;
            end
            if isfield(event,'end_ind')
                obj.end_ind = event.end_ind;
            end
            if isfield(event,'start_time')
                obj.start_time = event.start_time;
            end
            if isfield(event,'end_time')
                obj.end_time = event.end_time;
            end
            if isfield(event,'ended_manually')
                obj.ended_manually = event.ended_manually;
            end
            if isfield(event,'start_file')
                obj.start_file = event.start_file;
            end
            if isfield(event,'end_file')
                obj.end_file = event.end_file;
            end
            if isfield(event,'sequence')
                obj.sequence = event.sequence;
            end
            if isfield(event,'predicted_levels')
                obj.predicted_levels = event.predicted_levels;
            end
            if isfield(event,'predicted_levels_stdev')
                obj.predicted_levels_stdev = event.predicted_levels_stdev;
            end
            if isfield(event,'level_means')
                obj.level_means = event.level_means;
            end
            if isfield(event,'level_medians')
                obj.level_medians = event.level_medians;
            end
            if isfield(event,'level_stds')
                obj.level_stds = event.level_stds;
            end
            if isfield(event,'level_timing')
                obj.level_timing = event.level_timing;
            end
            if isfield(event,'level_finding_params')
                obj.level_finding_params = event.level_finding_params;
            end
            if isfield(event,'level_alignment')
                obj.level_alignment = event.level_alignment;
            end
            if isfield(event,'pulses')
                obj.pulses = event.pulses;
            end
            if isfield(event,'voltage')
                obj.voltage = event.voltage;
            end
            if isfield(event,'temp')
                obj.temp = event.temp;
            end
            if isfield(event,'KCl_molarity')
                obj.KCl_molarity = event.KCl_molarity;
            end
            if isfield(event,'ADPNP_molarity')
                obj.ADPNP_molarity = event.ADPNP_molarity;
            end
            if isfield(event,'open_pore_current')
                obj.open_pore_current = event.open_pore_current;
            end
        end
        
        function save(obj)
            % save data for a molecule object so it can be re-instantiated
            % later
            warning off MATLAB:structOnObject
            molecule_data = struct(obj);
            warning on MATLAB:structOnObject
            % autosaves to
            % /Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/__date__/__date__file__starttime_mol.mat
            file_dest = [obj.start_file(1:42) 'Analysis' obj.start_file(47:end-4) '_' num2str(round(obj.start_time)) 's_mol.mat'];
            save(file_dest, 'molecule_data');
            display(['Saved molecule data in ' file_dest])
        end
        
        function addData(obj, event)
            % add more to this molecule, probably from a different file
            % do a check to see if the end is near the end of a file... if
            % not, display a warning
            sigdata = SignalData(obj.start_file);
            if obj.end_time < sigdata.tend - 0.001 % not really at end...
                display('Warning: you are adding data to a molecule that appears to be complete.')
            end
            try
                obj.end_time = event.end_time;
                obj.end_ind = event.end_ind;
                obj.ended_manually = event.ended_manually;
                obj.end_file = event.end_file;
                if isempty(obj.voltage) || isnan(obj.voltage)
                    if ~isempty(event.voltage) && ~isnan(event.voltage)
                        obj.voltage = event.voltage;
                    end
                end
                if isempty(obj.open_pore_current) || isnan(obj.open_pore_current)
                    if ~isempty(event.open_pore_current) && ~isnan(event.open_pore_current)
                        obj.open_pore_current = event.open_pore_current;
                    end
                end
            catch
                display('Error: tried to add data that is incomplete.');
            end
        end
        
        % plot functions
        function f = plot_current(obj, varargin)
            % plot_current(filter, trange, sampling_frequency) optional inputs
            
            % handle inputs
            in = obj.parseOptionalInputs(varargin{:});
            filter_freq = in.filter;
            trange = in.trange;
            
            n = 5000; % max number of points to plot per signal
            d = util.doLoadMoleculeData(obj, n, 'minmax', filter_freq, trange);
            v = util.doLoadMoleculeViewData(obj, n, trange);
            
            % display as positive
            d(:,2) = abs(d(:,2));
            v(:,2) = abs(v(:,2));
            
            cmap = get(groot,'defaultAxesColorOrder');
            
            f = figure();
            
            % adjust the axes if a specific time window was requested,
            % since loading returns time starting from zero
            if ~isempty(trange)
                d(:,1) = d(:,1)+trange(1);
                v(:,1) = v(:,1)+trange(1);
            end
            
            % plot
            plot(v(:,1),v(:,2),'Color',0.8*ones(1,3)) % concatenates if off by a few points in length
            hold on
            plot(d(:,1),d(:,2),'Color','k')
            xlim([d(1,1), d(end,1)]);
            buff = 100;
            ylimits = [max(0,min(d(buff:end-buff,2))-20), max(d(buff:end-buff,2))+30];
            ylim(ylimits);
            % include pulses if we have them
            sigdata = SignalData(obj.end_file);
            if (sigdata.nsigs == 3)
                % check if we need a second file
                obj.pulses = obj.getPulseTiming(obj, 4, trange);
                line(ones(2,1)*obj.pulses,ylimits'*ones(size(obj.pulses)),'Color','r') % red vertical line at each pulse
            end
            % pretty it up
            title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) ...
                '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4)],'FontSize',24)
            if ~strcmp(obj.start_file, obj.end_file)
                title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) ...
                    '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ...
                    ' - ' obj.end_file(end-7:end-4)],'FontSize',24)
            end
            xlabel('Time (s)')
            ylabel('Current (pA)')
            annotation('textbox', [0.78 0.87 0 0], 'String', [num2str(filter_freq) 'Hz'], 'FontSize', 20, 'Color', 'k');
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_quick(obj, varargin)
            % plot_quick(trange) optional input
            
            % handle inputs
            in = obj.parseOptionalInputs(varargin{:});
            trange = in.trange;
            
            v = util.doLoadMoleculeViewDataQuick(obj, trange);
            
            % display as positive
            v(:,2) = abs(v(:,2));
            
            cmap = get(groot,'defaultAxesColorOrder');
            
            f = figure();
            
            % adjust the axes if a specific time window was requested,
            % since loading returns time starting from zero
            if ~isempty(trange)
                v(:,1) = v(:,1)+trange(1);
            end
            
            % plot
            plot(v(:,1),v(:,2),'Color',0.6*ones(1,3)) % concatenates if off by a few points in length
            xlim([v(1,1), v(end,1)]);
            buff = 100;
            ylimits = [max(0,min(v(buff:end-buff,2))-20), max(v(buff:end-buff,2))+30];
            ylim(ylimits);
            
            % pretty it up
            title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) ...
                '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4)],'FontSize',24)
            if ~strcmp(obj.start_file, obj.end_file)
                title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) ...
                    '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ...
                    ' - ' obj.end_file(end-7:end-4)],'FontSize',24)
            end
            xlabel('Time (s)')
            ylabel('Current (pA)')
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_squiggle(obj, varargin)
            % plot squiggle, optional arguments determine level refinement
            f = figure();
            if numel(varargin)>0
                l = obj.get_robust_levels(varargin{:});
            else
                l.level_means = obj.level_means;
                l.level_stds = obj.level_stds;
            end
            errorbar(1:numel(l.level_means),abs(l.level_means),l.level_stds,'o-')
            
            title(['Squiggle data, ' num2str(round(obj.voltage)) 'mV, ' num2str(obj.temp) '°C'],'FontSize',24)
            xlabel('Level')
            ylabel('Mean current (pA)')
            xlim([0, numel(l.level_means)+1])
            annotation('textbox', [0.8 0.87 0 0], 'String', ...
                [obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ...
                ' [' num2str(round(obj.start_time)) ',' num2str(round(obj.end_time)) '] ' ...
                ' p=' num2str(obj.level_finding_params.p)], ...
                'FontSize', 20);
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_levels(obj,varargin)
            % plot_levels(trange) optional input.  if trange is specified,
            % this will use the actual sampling for the plot.
            
            f = obj.plot_quick(varargin{:});
            
            hold on
            line(obj.level_timing'-obj.level_timing(1,1),(abs(obj.level_means)*[1,1])','LineWidth',2);
            robust_levs = obj.get_robust_levels(varargin{:});
            line(robust_levs.level_timing'-robust_levs.level_timing(1,1),(abs(robust_levs.level_means)*[1,1])','LineWidth',3,'Color','r');
            title(['Found levels, ' num2str(round(obj.voltage)) 'mV, ' num2str(obj.temp) '°C'],'FontSize',24)
            xlabel('Time (s)')
            ylabel('Current (pA)')
            annotation('textbox', [0.8 0.82 0 0], 'String', ...
                [obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ' p=' ...
                num2str(obj.level_finding_params.p) ' s=' num2str(obj.level_finding_params.sampling) ...
                'Hz'], 'FontSize', 20);
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_alignment(obj)
            % plot the alignment information
            % if alignment hasn't been done, do it
            if ~isfield(obj.level_alignment,'model_level_assignment')
                %obj.do_iterative_scaling_alignment;
                display('Alignment not completed!')
                return;
            end
            f = figure();
            
            mod_inds = obj.level_alignment.model_level_assignment;
            mod_type = obj.level_alignment.level_type; % 1 is normal, 2 is noise, 3 is deep block
            lvl_accum = obj.level_alignment.model_levels_measured_mean_currents; % mean of level currents for each level assigned to a given model level
            lvls = obj.level_means;
            
            subplot(3,1,1);
            plot(obj.predicted_levels,'o-','LineWidth',2)
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
            xx = 1:numel(obj.level_means);
            plot(xx(mod_type==2),abs(obj.level_means(mod_type==2)),'rx','MarkerSize',10)
            plot(xx(mod_type==3),abs(obj.level_means(mod_type==3)),'go','MarkerSize',10)
            ylabel('Current (pA)')
            xlabel('Measured level')
            xlim([0 numel(lvls)+1])
            set(gca,'FontSize',20)
            title('Matching each measured level to a model state')
            
            subplot(3,1,3);
            imagesc(1.03.^(obj.level_alignment.P/2)) % scales so image shows up well
            hold on
            plot(1:numel(lvls),obj.level_alignment.ks,'r','LineWidth',3);
            title('Probability Matrix')
            ylabel('State')
            xlabel('Level Number')
            set(gca,'FontSize',20)
            
            set(f,'position',[-1049, -382, 1050, 1585]);
        end
        
        function f = plot_enzyme_location(obj)
            
            f = figure();
            logic = obj.level_alignment.level_type==1;
            plot(obj.level_alignment.model_level_assignment(logic),'o-')
            xlabel('Measured level')
            xlim([0, numel(obj.level_means)+1])
            
            t = cumsum(diff(obj.level_timing,1,2));
            logic = obj.level_alignment.level_type==1;
            plot(t(logic), obj.level_alignment.model_level_assignment(logic),'o-')
            xlabel('Time (s)')
            xlim([0, t(end)+1])
            
            title(['Enzyme position, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) '°C'],'FontSize',24)
            ylabel(['Location on DNA strand' sprintf('\n') 'Model level number'])
            annotation('textbox', [0.12 0.87 0 0], 'String', ...
                [obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ...
                ' [' num2str(round(obj.start_time)) ',' num2str(round(obj.end_time)) '] ' ...
                num2str(obj.level_finding_params.filter) 'Hz, p=' num2str(obj.level_finding_params.p)], ...
                'FontSize', 20);
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
            
        end
        
        function p = plot_pulse_scatter(obj, varargin)
            % plot_pulse_scatter(filter, pA_threshhold)
            
            % handle inputs
            in = obj.parseOptionalInputs(varargin{:});
            filter = in.filter;
            pA_threshhold = in.threshhold;
            
            f = figure();
            hold on
            set(gca,'fontsize',20)
            ylabel('Current (pA)')
            xlabel('Time (ms)')
            title('Pulses overlaid: current scatter plot')
            set(f,'Position',[-965, 272, 852, 862])
            
            sigdata = SignalData(obj.start_file);
            filtname = sprintf('Low-pass (%d Hz)', filter);
            fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
            
            fprintf(['Pulses with a minimum ' num2str(pA_threshhold) 'pA change: '])
            p = [];
            for i = 1:numel(obj.pulses)
                % grab data
                ms_before = 1;
                ms_after = 10;
                d = sigdata.get([obj.pulses(i)-ms_before/1000, obj.pulses(i)+ms_after/1000]/sigdata.si);
                %time = d(:,1)-m.pulses(i);
                time = linspace(-1*ms_before, ms_after, size(d,1));
                current = d(:,fsigs(1))*1000;
                % subtract starting level
                I_before = mean(current(round((-0.0003/sigdata.si:-0.0002/sigdata.si)+ms_before/1000/sigdata.si)));
                I_after = mean(current(round((0.0012/sigdata.si:0.0025/sigdata.si)+ms_before/1000/sigdata.si)));
                if abs(I_after-I_before)>pA_threshhold
                    plot(time+randn(1)*0.002,current-I_before,'k.','MarkerSize',0.5)
                    drawnow;
                    fprintf([num2str(i) ', '])
                    p = [p, i];
                end
            end
            
            fprintf('\b\b.\n')
            display(['Number = ' num2str(numel(p)) ', out of a total of ' num2str(numel(obj.pulses)) ' in the entire molecule.'])
            
        end
        
        function f = plot_pulse(obj, varargin)
            % plot_pulse_scatter(pulse_number)
            
            % handle inputs
            in = obj.parseOptionalInputs(varargin{:});
            filter = in.filter;
            num = in.pulsenum;
            
            f = figure(1);
            clf(1)
            hold on
            set(gca,'fontsize',20)
            ylabel('Current (pA)')
            xlabel('Time (ms)')
            title([obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ': pulses ' num2str(num)])
            set(f,'Position',[-965, 272, 852, 862])
            
            sigdata = SignalData(obj.start_file);
            filtname = sprintf('Low-pass (%d Hz)', filter);
            fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
            
            twindow = [-0.05, 0.05];
            if numel(num)==1
                d = sigdata.getByTime(obj.pulses(num)+twindow);
                plot((d(:,1)-obj.pulses(num))*1000,d(:,fsigs(1))*1000,'k')
                plot(0,mean(d(1:100,fsigs(1)))*1000+50,'r*')
                ylim([0, 500])
            else
                for n = num
                    d = sigdata.getByTime(obj.pulses(n)+twindow);
                    plot((d(:,1)-obj.pulses(n))*1000,(d(:,fsigs(1))-mean(d(1:size(d,1)/3,fsigs(1))))*1000)
                end
                xlim(twindow*1000)
            end
        end
        
        function levels = do_level_analysis(obj, varargin)
            
            % handle inputs
            in = obj.parseOptionalInputs(varargin{:});
            filter = in.filter;
            downsample_frequency = in.sample;
            p = in.plevels;
            
            % do level finding
            levels = find_levels(obj, filter, downsample_frequency, p);
            %levels = find_levels_ks(obj, filter, downsample_frequency, p);
            
            % package the data
            obj.level_means = cellfun(@(x) x.current_mean, levels);
            obj.level_medians = cellfun(@(x) x.current_median, levels);
            obj.level_stds = cellfun(@(x) x.current_std, levels);
            obj.level_timing = [cellfun(@(x) x.start_time, levels), cellfun(@(x) x.end_time, levels)];
            obj.level_finding_params.p = p;
            obj.level_finding_params.sampling = downsample_frequency;
            obj.level_finding_params.filter = filter;
            
        end
        
        function do_iterative_scaling_alignment(obj)
            % iterates between model scaling and level alignment until
            % convergence is reached.  helps improve model scaling.
            
            % get the initial predicted levels from oxford
            levs = abs(obj.level_means);
            hicut = 0.57 * abs(obj.open_pore_current);
            lowcut = 0.1 * abs(obj.open_pore_current);
%             [model_levels, model_levels_std] = ...
%                 get_model_levels_oxford(obj.sequence, levs(levs>lowcut & levs<hicut), abs(obj.open_pore_current), abs(obj.voltage), obj.temp);
            [model_levels, model_levels_std,~,~] = ...
                get_model_levels_M2(obj.sequence, levs(levs>lowcut & levs<hicut));
            
            % save the initial scaling
            obj.predicted_levels = model_levels';
            obj.predicted_levels_stdev = model_levels_std';
            
            % iterate alignment and scaling until convergence is achieved
            stop_criterion = 0.05;
            lsqdist = [];
            dfit = 1;
            while dfit > stop_criterion && numel(lsqdist) < 20
                % do a level alignment
                obj.do_level_alignment;
                % calculate least-squares distance per measured level
                logic = ~isnan(obj.level_alignment.model_levels_measured_mean_currents); % these levels are not missing
                lsqdist(end+1) = sum((obj.level_alignment.model_levels_measured_mean_currents(logic) - obj.predicted_levels(logic)).^2) / sum(logic);
                if numel(lsqdist)>1
                    dfit = abs(lsqdist(end)-lsqdist(end-1));
                end
                if sum(logic)<2
                    display('Cannot do iterative fitting, too far off.')
                    break;
                end
                % do a least-squares fit
                f = fit(obj.level_alignment.model_levels_measured_mean_currents(logic),obj.predicted_levels(logic), ...
                        'poly1','Weights',obj.level_alignment.model_levels_measured_total_duration(logic),'Robust','bisquare');
                % invert to get the scale and offset corrections and update the
                % predicted levels
                obj.predicted_levels = (obj.predicted_levels-f.p2)/f.p1;
            end
            
            if numel(lsqdist)==20
                display('Iterative scaling alignment warning: convergence not reached in 20 iterations.')
            end
            
        end
        
        function do_level_alignment(obj)
            % align the levels to the predicted model levels
            % error if we are missing something
            if isnan(obj.level_means)
                display('Error: could not align levels.  No levels extracted from the data!');
                return;
            end
            if isnan(obj.predicted_levels)
                display('Error: could not align levels.  No predicted sequence levels!');
                return;
            end
            obj.level_alignment = struct(); % clear any previous alignment
            % get reasonable levels
            [mod_inds, mod_type, lvl_accum, P, ks] = align_fb(obj.predicted_levels, ...
                obj.predicted_levels_stdev, abs(obj.level_means), diff(obj.level_timing,1,2), 0.18*abs(obj.open_pore_current));
            obj.level_alignment.model_level_assignment = mod_inds;
            obj.level_alignment.level_type = mod_type; % 1 is normal, 2 is noise, 3 is deep block
            obj.level_alignment.model_levels_measured_mean_currents = lvl_accum; % mean of level currents for each level assigned to a given model level
            obj.level_alignment.model_levels_measured_total_duration = accumarray(mod_inds(mod_type~=2), ...
                diff(obj.level_timing(mod_type~=2,:),1,2),size(obj.predicted_levels),@sum,nan); % total time in each level
            obj.level_alignment.P = P; % probabilities in the alignment matrix
            obj.level_alignment.ks = ks; % state index in the alignment matrix
            
            % also store the level means and timing after alignment, by
            % which i mean: combine adjacent "stay" levels
            lmean = [];
            lmed = [];
            lstd = [];
            ltime = [];
            level_means_no_noise = abs(obj.level_means(mod_type==1));
            level_medians_no_noise = abs(obj.level_medians(mod_type==1));
            level_stds_no_noise = obj.level_stds(mod_type==1);
            level_timing_no_noise = obj.level_timing(mod_type==1,:);
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
            obj.level_alignment.level_means = lmean';
            obj.level_alignment.level_medians = lmed';
            obj.level_alignment.level_stds = lstd';
            obj.level_alignment.level_timing = ltime - ltime(1,1); % zero at beginning of molecule
            
        end
        
        function update_pulse_timing(obj)
            % update the pulse timings
            sigdata = SignalData(obj.start_file);
            if (sigdata.nsigs == 3)
                p = obj.getPulseTiming();
                obj.pulses = p;
            end
        end
        
        function levs = get_robust_levels(obj, varargin)
            % combines levels it thinks are stays
            % eliminates levels it thinks are noise
            % but all without actually aligning to a model set of levels
            
%             % use Kolmogorov-Smirnov to combine stays
%             % empirically determined level is the KS parameter < 0.07
%             % means the level is the same
%             
            % handle inputs
            in = obj.parseOptionalInputs(varargin{:});
            p_s = in.pstay;
            p_n = in.pnoise;
%             
%             % grab the levels
%             levs.m = obj.level_means;
%             levs.s = obj.level_stds;
%             levs.t = obj.level_timing;
%             
%             % note the obvious noise levels
%             noise = obj.level_stds>8 | abs(obj.level_means)<abs(obj.open_pore_current)*0.08;
%             
%             % initial setting of first level
%             if ~noise(1)
%                 goodlevs.m(1) = levs.m(1);
%                 goodlevs.s(1) = levs.s(1);
%                 goodlevs.t(1,:) = levs.t(1,:);
%             else
%                 goodlevs.m(1) = levs.m(2);
%                 goodlevs.s(1) = levs.s(2);
%                 goodlevs.t(1,1) = levs.t(1,1);
%                 goodlevs.t(1,2) = levs.t(2,2);
%             end
%             
%             % go through levels
%             for i = 2:numel(levs.m)
%                 
%                 % check if the level is noise
%                 if noise(i)
%                     % ignore it and keep the stats unchanged
%                     goodlevs.t(end,2) = levs.t(i,2);
%                 else
%                     % check for a possible stay by looking at means and std
%                     if abs(levs.m(i)-goodlevs.m(end)) < sqrt(goodlevs.s(end)^2+levs.s(i)^2)
%                         % check if the level is a stay using KS test
%                         window = [goodlevs.t(end,1), levs.t(i,2)];
%                         % display(['loading ' num2str(diff(window)) ' sec'])
%                         tic;
%                         data = util.doLoadMoleculeData(obj, diff(window)*5000, 'none', 2000, window); % 2 kHz filter, minmax
%                         toc;
%                         ind = round(diff(goodlevs.t(end,:))/diff(window)*size(data,1));
%                         lastlev = randsample(data(ind:end,2),min(1000,size(data,1)-ind));
%                         thislev = randsample(data(1:ind,2),min(1000,ind));
%                         [~,~,k] = kstest2(lastlev,thislev);
%                         %display(num2str(k))
%                         if k < 0.2%0.07
%                             % this is the same level
%                             goodlevs.m(end) = (levs.m(i)*diff(levs.t(i,:))+goodlevs.m(end)*diff(goodlevs.t(end,:)))/(diff(window));
%                             goodlevs.s(end) = (levs.s(i)*diff(levs.t(i,:))+goodlevs.s(end)*diff(goodlevs.t(end,:)))/(diff(window));
%                             goodlevs.t(end,2) = levs.t(i,2); % update end time for combined level
%                             % levdata = util.doLoadMoleculeData(obj, 100, 'pointwise', 2000, goodlevs.t(end,:)); % 100pt pointwise downsample
%                             % goodlevs.m(end) = nanmean(levdata(:,2)); % update mean
%                             % goodlevs.s(end) = nanstd(levdata(:,2)); % update standard deviation
%                             
%                         else
%                             % this is a new level
%                             goodlevs.m(end+1) = levs.m(i);
%                             goodlevs.s(end+1) = levs.s(i);
%                             goodlevs.t(end+1,:) = levs.t(i,:);
%                         end
%                     else
%                         % this is a new level
%                         goodlevs.m(end+1) = levs.m(i);
%                         goodlevs.s(end+1) = levs.s(i);
%                         goodlevs.t(end+1,:) = levs.t(i,:);
%                     end
%                 end
%                 
%             end
%             
%             clear levs
%             levs.level_means = goodlevs.m'; % replace initial molecule levels with these new ones
%             levs.level_stds = goodlevs.s';
%             levs.level_timing = goodlevs.t;
            
            lev = abs(obj.level_means);
            stdv = abs(obj.level_stds);
            timing = obj.level_timing;
            dur = obj.level_timing(:,2)-obj.level_timing(:,1);
            % remove levels with huge standard deviation and deep blocks
            lev = lev(obj.level_stds<10 & abs(obj.level_means)>abs(obj.open_pore_current)*0.08);
            stdv = stdv(obj.level_stds<10 & abs(obj.level_means)>abs(obj.open_pore_current)*0.08);
            timing = timing(obj.level_stds<10 & abs(obj.level_means)>abs(obj.open_pore_current)*0.08,:);
            dur = dur(obj.level_stds<10 & abs(obj.level_means)>abs(obj.open_pore_current)*0.08);
            tau = 1e-3;
            levs.level_means(1) = lev(1);
            levs.level_stds(1) = stdv(1);
            levs.level_timing(1,:) = obj.level_timing(1,:);
            
            for i = 2:numel(lev)
                % probability this is a stay
                stdev = levs.level_stds(end);
                div = 0.00005;
                p_stay = mean([ normpdf(lev(i),lev(i-1)+div,stdev/sqrt(timing(i)/(100*tau))), normpdf(lev(i),lev(i-1)-div,stdev/sqrt(timing(i)/(100*tau))) ]) * 2*div ... % prob of distribution overlap
                    * 2^(-max(0,abs(stdev-stdv(i))-0.1)/(0.1*min(stdv(i),stdev))); % factor to make stays unlikely for different standard deviations
                % probability this is noise
                if i<numel(lev)
                    p_noise = mean([ exppdf(dur(i)+div,tau), exppdf(dur(i)-div,tau) ]) * 2*div ... % short enough
                         * (1 - p_stay) ... % it's not the same level
                         * mean([ normpdf(lev(i+1),lev(i-1)+div,stdev), normpdf(lev(i+1),lev(i-1)-div,stdev) ]) * 2*div; % it goes back to the same level
                else
                    p_noise = 0;
                end
                % do stuff
                %line(levs.level_timing',(abs(levs.level_means)*[1,1])','LineWidth',2,'Color','r');
                %line(obj.level_timing(i,:)',(abs(obj.level_means(i))*[1,1])','LineWidth',2,'Color','c');
                %display(['p_noise = ' num2str(log10(p_noise)) ', p_stay = ' num2str(log10(p_stay))])
                %pause();
                if log10(p_noise) > p_n
                    levs.level_timing(end,2) = timing(i,2); % update end of level
                elseif log10(p_stay) > p_s
                    levs.level_means(end) = ( levs.level_means(end)*diff(levs.level_timing(end,:)) + lev(i)*dur(i) ) ...
                        / (diff(levs.level_timing(end,:))+dur(i)); % mean over the whole level
                    levs.level_stds(end) = ( levs.level_stds(end)*diff(levs.level_timing(end,:)) + stdv(i)*dur(i) ) ...
                        / (diff(levs.level_timing(end,:))+dur(i)); % std over the whole level
                    levs.level_timing(end,2) = timing(i,2); % update end of level
                else
                    levs.level_means(end+1,1) = lev(i);
                    levs.level_stds(end+1,1) = stdv(i);
                    levs.level_timing(end+1,:) = timing(i,:);
                end
                p1(i) = p_stay;
                p2(i) = p_noise;
            end
        end
        
        function n = get_alignment_stats(obj, model_levels_of_interest)
            % return alignment statistics pertaining to coverage and
            % forward and backward steps
            
            if nargin < 2
                model_levels_of_interest = 1:numel(obj.predicted_levels);
            end
            if strcmp(model_levels_of_interest,'end')
                model_levels_of_interest(2) = numel(obj.predicted_levels);
            end
            
            % go through the levels and count
            assignments = obj.level_alignment.model_level_assignment;
            real_levels = obj.level_alignment.level_type == 1 & ismember(assignments, model_levels_of_interest);
            ll = assignments(real_levels);
            
            n.skips = 0;
            n.stays = 0;
            n.forward = 0;
            n.back = 0;
            
            for i = 2:numel(ll)
                delta = ll(i)-ll(i-1);
                switch delta
                    case -1
                        n.back = n.back + 1;
                    case 0
                        n.stays = n.stays + 1;
                    case 1
                        n.forward = n.forward + 1;
                    case 2
                        n.skips = n.skips + 1;
                    otherwise
                        if ll(i)-ll(i-1) < -1
                            n.back = n.back + abs(ll(i)-ll(i-1));
                        elseif ll(i)-ll(i-1) > 2
                            n.skips = n.skips + abs(ll(i)-ll(i-1)) - 1;
                        end
                end
            end
            
            n.noise = sum(obj.level_alignment.level_type == 2);
            n.deepblock = sum(obj.level_alignment.level_type == 3);
            n.model_levels_covered = sum(~isnan(obj.level_alignment.model_levels_measured_mean_currents));
            n.model_levels_never_covered = sum(isnan(obj.level_alignment.model_levels_measured_mean_currents(1:find(~isnan(obj.level_alignment.model_levels_measured_mean_currents),1,'last'))));
        end
        
        function kappa = pulse_correlation_metric(obj)
            % a metric for correlation between pulses and level changes
            
            if isnan(obj.pulses)
                kappa = nan;
                return;
            end
            if isempty(obj.pulses)
                kappa = nan;
                return;
            end
            
            % cumulative distribution function for level durations
            n = obj.get_robust_levels(-6,-5);
            tau = mean(diff(n.level_timing,1,2));
            c_d_f = @(t) 1 - exp(-1*t/tau);
            
            % get pulse timings and level change timings
            level_change_timings = n.level_timing(:,2);
            pulse_timings = obj.pulses;
            
            % only look at level changes during pulsing time
            i1 = find(level_change_timings>pulse_timings(1),1,'first');
            i2 = find(level_change_timings>pulse_timings(end),1,'first');
            level_change_timings = level_change_timings(i1:i2);
            
            % calculate p value for each pulse
            p = nan(1,numel(pulse_timings));
            for i = 1:numel(pulse_timings)-1
                ind = find(level_change_timings>pulse_timings(i)&level_change_timings<pulse_timings(i+1),1,'first'); % closest level change
                if ~isempty(ind)
                    p(i) = c_d_f(level_change_timings(ind)-pulse_timings(i));
                end
            end
            
            % calculate kappa
            bins = ceil(numel(pulse_timings)/10);
            bincenters = 1/bins/2 : 1/bins : 1 - 1/bins/2;
            n = hist(p, bincenters);
            kappa = sqrt( std(n)^2*bins / numel(p)^2 );
            
            % plot if no arguments passed out
            if nargout==0
                figure(1)
                clf(1)
                bar(bincenters,n)
                xlabel('p value')
                ylabel('Number of pulses')
                title(['kappa = ' num2str(kappa,2)])
                set(gca,'fontsize',20)
                display(kappa)
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function p = getPulseTiming(obj, varargin) % NEEDS UPDATE!
            % getPulseTiming(sigdata, channel, trange)
            % check if we've done this already
            if ~isnan(obj.pulses)
                p = obj.pulses;
                return;
            end
            % inputs
            if numel(varargin)==0
                if ~strcmp(obj.start_file, obj.end_file)
                    sigdata = [SignalData(obj.start_file), SignalData(obj.end_file)];
                    trange = [obj.start_time, sigdata.tend; 0, obj.end_time];
                else
                    sigdata = SignalData(obj.start_file);
                    trange = [obj.start_time, obj.end_time];
                end
                channel = 4;
            else
                sigdata = varargin{1};
                channel = varargin{2};
                trange = varargin{3};
            end
            % find possible pulses quickly using a derivative
            for j = 1:numel(sigdata)
                d = sigdata(j).getViewData(trange(j,:));
                dt = d(2,1)-d(1,1);
                pulsedata = d(:,channel);
                num = ceil(0.001/dt);
                difference = diff(medfilt1(pulsedata,num));
                [~,candidates] = findpeaks(difference.*(difference>0),'MinPeakHeight',1,'MinPeakDist',50);
                candidates = trange(j,1) + candidates*dt; % convert from index to time
                % refine pulse timing
                p = [];
                for i = 1:numel(candidates)
                    ind = (candidates(i)-2e-3)/sigdata.si;
                    p(end+1) = sigdata.findNext(@(x) x(:,4)>0.1, ind) * sigdata.si;
                end
            end
        end
        
        function in = parseOptionalInputs(obj, varargin)
            % use inputParser to figure out what all these inputs mean
            p = inputParser;
            
            % defaults and checks
            defaultFilterFreq = 1000;
            defaultSampleFreq = 5000;
            checkFilterFreq = @(x) all([isnumeric(x),x>=10,x<=20000]);
            
            defaultTimeRange = [];
            checkTimeRange = @(x) all([x(1)>=0, x(1)<x(2), diff(x)<=obj.level_timing(end,2)-obj.level_timing(1,1)]);
            
            defaultPstay = -4;
            defaultPnoise = -10;
            defaultPlevels = -15;
            checkP = @(x) all([isnumeric(x), x<=0]);
            
            % set up the inputs
            addOptional(p,'filter',defaultFilterFreq,checkFilterFreq);
            addOptional(p,'sample',defaultSampleFreq,checkFilterFreq);
            addOptional(p,'trange',defaultTimeRange,checkTimeRange);
            addOptional(p,'pstay',defaultPstay,checkP);
            addOptional(p,'pnoise',defaultPnoise,checkP);
            addOptional(p,'plevels',defaultPlevels,checkP);
            addOptional(p,'pulsenum',[],@isnumeric);
            addOptional(p,'threshhold',5,@isnumeric);
            
            % parse
            parse(p,varargin{:});
            in = p.Results;
        end
        
    end
    
end