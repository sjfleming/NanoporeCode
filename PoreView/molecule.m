classdef molecule < handle & matlab.mixin.SetGet
    %MOLECULE object to aid in analysis of nanopore sequencing data
    %   ...
    
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
            catch
                display('Error: tried to add data that is incomplete.');
            end
        end
        
        % plot functions
        function f = plot_current(obj, varargin)
            % plot_current(filter, trange, sampling_frequency) optional inputs
            n = 10000; % max number of points to plot per signal
            if numel(varargin) == 0
                filter_freq = 1000;
                trange = [obj.start_time, obj.end_time];
            elseif numel(varargin) == 1
                filter_freq = varargin{1};
                trange = [obj.start_time, obj.end_time];
            elseif numel(varargin) == 2
                filter_freq = varargin{1};
                trange = varargin{2};
                if ~(trange(1) >= obj.start_time && trange(2) <= obj.end_time) && ~strcmp(obj.start_file, obj.end_file)
                    display('Requested time range is outside the range of the molecule.');
                    return;
                end
            elseif numel(varargin) == 3
                filter_freq = varargin{1};
                trange = varargin{2};
                if ~(trange(1) >= obj.start_time && trange(2) <= obj.end_time) && ~strcmp(obj.start_file, obj.end_file)
                    display('Requested time range is outside the range of the molecule.');
                    return;
                end
                n = diff(trange) * varargin{3};
            end
            cmap = get(groot,'defaultAxesColorOrder');
            
            f = figure();
            sigdata = SignalData(obj.start_file);
            filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter_freq);
            fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter_freq),filtname);
            
            % check to see if we are using one or two files
            if (strcmp(obj.start_file, obj.end_file) || trange(2) > trange(1))
                current = util.downsample_pointwise(sigdata,fsigs(1),trange,n)*1000; % filtered data, in pA
                current_raw = util.downsample_minmax(sigdata,2,trange,n)*1000; % in pA
                time = linspace(trange(1), trange(2), numel(current));
            else
                sigdata2 = SignalData(obj.end_file);
                n1 = (sigdata.tend - trange(1)) / (trange(2) + sigdata.tend - trange(1)) * n;
                n2 = (trange(2)) / (trange(2) + sigdata.tend - trange(1)) * n;
                trange1 = [trange(1), sigdata.tend];
                current1 = util.downsample_pointwise(sigdata,fsigs(1),trange1,n1)*1000; % filtered data, in pA
                current_raw1 = util.downsample_minmax(sigdata,2,trange1,n1)*1000; % in pA
                trange2 = [0, trange(2)];
                fsigs2 = sigdata2.addVirtualSignal(@(d) filt_lpb(d,4,filter_freq),filtname);
                current2 = util.downsample_pointwise(sigdata2,fsigs2(1),trange2,n2)*1000; % filtered data, in pA
                current_raw2 = util.downsample_minmax(sigdata2,2,trange2,n2)*1000; % in pA
                current = [current1, current2];
                current_raw = [current_raw1, current_raw2];
                time = linspace(trange(1), sigdata.tend + trange(2), numel(current));
            end
            
            % plot
            plot(time,current_raw(1:numel(time)),'Color',0.8*ones(1,3)) % concatenates if off by a few points in length
            hold on
            plot(time,current,'Color','k')
            xlim([time(1), time(end)]);
            buff = 100;
            ylimits = [max(0,min(current_raw(buff:end-buff))-20), max(current_raw(buff:end-buff))+20];
            ylim(ylimits);
            if (sigdata.nsigs == 3)
                % check if we need a second file
                if (strcmp(obj.start_file, obj.end_file) || trange(2) > trange(1))
                    obj.pulses = obj.getPulseTiming(sigdata, 4, trange);
                else
                    trange1 = [trange(1), sigdata.tend];
                    trange2 = [0, trange(2)];
                    pulses1 = obj.getPulseTiming(sigdata, 4, trange1);
                    pulses2 = obj.getPulseTiming(sigdata2, 4, trange2);
                    obj.pulses = [pulses1, pulses2]; % concatenate
                end
                line(ones(2,1)*obj.pulses,ylimits'*ones(size(obj.pulses)),'Color','r') % red vertical line at each pulse
            end
            % pretty it up
            title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) ...
                '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4)],'FontSize',24)
            if (trange(2) < trange(1) && ~strcmp(obj.start_file, obj.end_file) )
                title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) ...
                    '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4) ...
                    ' - ' obj.end_file(end-7:end-4)],'FontSize',24)
            end
            xlabel('Time (s)')
            ylabel('Current (pA)')
            annotation('textbox', [0.8 0.87 0 0], 'String', [num2str(filter_freq) 'Hz'], 'FontSize', 20, 'Color', 'k');
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_squiggle(obj)
            f = figure();
            plot(1:numel(obj.level_means),obj.level_means,'o-')
            title(['Squiggle data, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) '°C'],'FontSize',24)
            xlabel('Level')
            ylabel('Mean current (pA)')
            xlim([0, numel(obj.level_means)+1])
            annotation('textbox', [0.8 0.87 0 0], 'String', ...
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
        
        function f = plot_levels(obj,varargin)
            % plot_levels(trange) optional input.  if trange is specified,
            % this will use the actual sampling for the plot.
            if (numel(varargin)==1)
                f = obj.plot_current(obj.level_finding_params.filter,varargin{1},obj.level_finding_params.sampling);
            else
                f = obj.plot_current(obj.level_finding_params.filter);
            end
            hold on
            line(obj.level_timing',(obj.level_means*[1,1])','LineWidth',3);
            title(['Found levels, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) '°C'],'FontSize',24)
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
                obj.do_level_alignment;
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
            plot(lvls)
            for i=1:numel(lvls)
                text(i,lvls(i),num2str(mod_inds(i)),'FontSize',14);
            end
            hold on
            xx = 1:numel(obj.level_means);
            plot(xx(mod_type==2),obj.level_means(mod_type==2),'rx','MarkerSize',10)
            plot(xx(mod_type==3),obj.level_means(mod_type==3),'go','MarkerSize',10)
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
        
        function p = plot_pulse_scatter(obj, varargin)
            % plot_pulse_scatter(filter, pA_threshhold)
            
            if numel(varargin)==0
                filter = 10000;
                pA_threshhold = 5;
            elseif numel(varargin)==1
                filter = varargin{1};
                pA_threshhold = 5;
            else
                filter = varargin{1};
                pA_threshhold = varargin{2};
            end
            
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
            
            filter = 10000;
            if numel(varargin)==0
                num = input('Input pulse number');
            else
                num = varargin{1};
            end
            
            f = figure();
            hold on
            set(gca,'fontsize',20)
            ylabel('Current (pA)')
            xlabel('Time (ms)')
            title([obj.start_file(end-27:end-20) ': pulses ' num2str(num)])
            set(f,'Position',[-965, 272, 852, 862])
            
            sigdata = SignalData(obj.start_file);
            filtname = sprintf('Low-pass (%d Hz)', filter);
            fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
            
            twindow = [-0.05, 0.05];
            if numel(num)==1
                d = sigdata.getByTime(obj.pulses(num)+twindow);
                plot(d(:,1),d(:,fsigs(1))*1000,'k')
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
        
        function levels = do_level_analysis(obj, filter, downsample_frequency, p)
            
            % access the data and filter and downsample
            sigdata = SignalData(obj.start_file);
            filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter);
            fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
            
            % if there are pulses, remove a time range around each to eliminate spikes
            % this is very necessary.  for some reason, spikes thow it off.
            if sigdata.nsigs>2
                obj.pulses = obj.getPulseTiming();
                fsigs = sigdata.addVirtualSignal(@(d) filt_rmrange(d,[obj.pulses'-1e-4, obj.pulses'+1e-3, nan(size(obj.pulses'))]),'spikes removed');
            end
            
            % check to see if we are using one or two files
            if strcmp(obj.start_file, obj.end_file)
                n = (obj.end_ind - obj.start_ind) * downsample_frequency*sigdata.si;
                trange = [obj.start_time, obj.end_time];
                current = util.downsample_pointwise(sigdata,fsigs(1),trange,n)*1000; % filtered data, in pA
                time = linspace(trange(1), trange(2), numel(current));
            else
                sigdata2 = SignalData(obj.end_file);
                n1 = (sigdata.ndata - obj.start_ind) * downsample_frequency*sigdata.si;
                n2 = obj.end_ind * downsample_frequency*sigdata2.si;
                trange1 = [obj.start_time, sigdata.tend];
                current1 = util.downsample_pointwise(sigdata,fsigs(1),trange1,n1)*1000; % filtered data, in pA
                trange2 = [0, obj.end_time];
                fsigs2 = sigdata2.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
                current2 = util.downsample_pointwise(sigdata2,fsigs2(1),trange2,n2)*1000; % filtered data, in pA
                current = [current1, current2];
                time = linspace(obj.start_time, sigdata.tend + obj.end_time, numel(current));
            end
            obj.voltage = round(mean(sigdata.get(obj.start_ind:obj.start_ind+100,3))); % voltage in mV to nearest mV
            
            % use Laszlo's level-finding algorithm to find levels
            pad = 10;
            levels = laszlo_levels([time(pad:end-pad)', current(pad:end-pad)'],p);
            
            % package the data
            obj.level_means = cellfun(@(x) x.current_mean, levels);
            obj.level_medians = cellfun(@(x) x.current_median, levels);
            obj.level_stds = cellfun(@(x) x.current_std, levels);
            obj.level_timing = [cellfun(@(x) x.start_time, levels), cellfun(@(x) x.end_time, levels)];
            obj.level_finding_params.p = p;
            obj.level_finding_params.filter = filter;
            obj.level_finding_params.sampling = downsample_frequency;
            
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
            [mod_inds, mod_type, lvl_accum, P, ks] = align_fb(obj.predicted_levels, obj.level_means, diff(obj.level_timing,1,2), 0.18*obj.open_pore_current);
            obj.level_alignment.model_level_assignment = mod_inds;
            obj.level_alignment.level_type = mod_type; % 1 is normal, 2 is noise, 3 is deep block
            obj.level_alignment.model_levels_measured_mean_currents = lvl_accum; % mean of level currents for each level assigned to a given model level
            obj.level_alignment.model_levels_measured_total_duration = accumarray(mod_inds(mod_type~=2),diff(obj.level_timing(mod_type~=2,:),1,2),size(obj.predicted_levels),@sum,nan); % total time in each level
            obj.level_alignment.P = P; % probabilities in the alignment matrix
            obj.level_alignment.ks = ks; % state index in the alignment matrix
        end
        
        function update_pulse_timing(obj)
            % update the pulse timings
            sigdata = SignalData(obj.start_file);
            if (sigdata.nsigs == 3)
                p = obj.getPulseTiming();
                obj.pulses = p;
            end
        end
        
    end
    
    methods (Access = private)
        
        function p = getPulseTiming(obj, varargin)
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
        
    end
    
end

