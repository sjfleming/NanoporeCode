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
        levels = nan;
        level_timing = nan;
        level_assignments = nan;
        pulses = nan;
        voltage = nan;
        temp = nan;
        KCl_molarity = nan;
        ADPNP_molarity = nan;
    end
    
    methods
        
        % constructor
        function obj = molecule(event)
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
            if isfield(event,'levels')
                obj.levels = event.levels;
            end
            if isfield(event,'level_timing')
                obj.level_timing = event.level_timing;
            end
            if isfield(event,'level_assignments')
                obj.level_assignments = event.level_assignments;
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
                obj.ended_manually = event.ended_manually;
                obj.end_file = event.end_file;
            catch
                display('Error: tried to add data that is incomplete.');
            end
        end
        
        % plot functions
        function f = plot_current(obj, varargin)
            % plot_current()
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
            end
            cmap = get(groot,'defaultAxesColorOrder');
            n = 1000; % max number of points to plot per signal
            f = figure();
            sigdata = SignalData(obj.start_file);
            filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter_freq);
            fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter_freq),filtname);
            current_raw = util.downsample(sigdata,2,trange,n)*1000; % in pA
            current = util.downsample(sigdata,fsigs(1),trange,n)*1000; % filtered data, in pA
            time = linspace(trange(1), trange(2), numel(current));
            obj.voltage = round(mean(sigdata.get(obj.start_ind:obj.start_ind+100,3))); % voltage in mV to nearest mV
            % check to see if we need a second file loaded
            if (trange(2) > sigdata.tend && ~strcmp(obj.start_file, obj.end_file) ) % if it overruns first file, and there are two files
                % concatenate from file 2
                sigdata2 = SignalData(obj.end_file);
                filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter_freq);
                fsigs = sigdata2.addVirtualSignal(@(d) filt_lpb(d,4,filter_freq),filtname);
                current_raw = [current_raw, util.downsample(sigdata2,2,trange,n)*1000]; % in pA
                current = [current, util.downsample(sigdata2,fsigs(1),trange,n)*1000]; % filtered data, in pA
                time = linspace(trange(1), trange(2), numel(current));
            end
            % plot
            plot(time,current_raw,'Color',0.8*ones(1,3))
            hold on
            plot(time,current,'Color',cmap(1,:))
            xlim(trange);
            buff = 100;
            ylimits = [max(0,min(current_raw(buff:end-buff))-20), max(current_raw(buff:end-buff))+20];
            ylim(ylimits);
            if (sigdata.nsigs == 3)
                obj.getPulseTiming(obj, sigdata, 4, trange); % updates obj.pulses variable
                line(ones(2,1)*obj.pulses,ylimits'*ones(size(obj.pulses)),'Color','r') % red vertical line at each pulse
            end
            % pretty it up
            title(['Full molecule sliding, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4)],'FontSize',24)
            xlabel('Time (s)')
            ylabel('Current (pA)')
            annotation('textbox', [0.85 0.87 0 0], 'String', [num2str(filter_freq) 'Hz'], 'FontSize', 20, 'Color', cmap(1,:));
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_squiggle(obj)
            f = figure();
            plot(1:numel(obj.levels),obj.levels,'o-')
            title(['Squiggle data, ' num2str(obj.voltage) 'mV, ' num2str(obj.temp) '°C : ' obj.start_file(end-27:end-20) '\_' obj.start_file(end-7:end-4)],'FontSize',24)
            xlabel('Level')
            ylabel('Mean current (pA)')
            annotation('textbox', [0.85 0.87 0 0], 'String', [num2str(filter_freq) 'Hz'], 'FontSize', 20, 'Color', cmap(1,:));
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
            set(gcf,'Position',[-1048, 803, 1000, 400]) % size the figure
            set(gcf,'Color',[1 1 1])
            set(gca,'FontSize',24)
        end
        
        function f = plot_alignment_matrix(obj)
            f = figure();
            
        end
        
    end
    
    methods (Access = private)
        
        function getPulseTiming(obj, sigdata, channel, trange)
            % check if we've done this already
            if ~isnan(obj.pulses)
                return;
            end
            % find possible pulses quickly using a derivative
            pds = 1000;
            d = sigdata.getViewData(trange);
            pulsedata = d(:,channel);
            difference = diff(pulsedata);
            [~,candidates] = findpeaks(difference.*(difference>0),'MinPeakHeight',2,'MinPeakDist',50);
            candidates = trange(1) + candidates/pds; % convert from index to time
            % refine pulse timing
            p = [];
            for i = 1:numel(candidates)
                ind = (candidates(i)-10e-3)/sigdata.si;
                p(end+1) = sigdata.findNext(@(x) x(:,4)>1, ind);
                p(end) = p(end) * sigdata.si; % convert from index to time
            end
            obj.pulses = p;
        end
        
    end
    
end

