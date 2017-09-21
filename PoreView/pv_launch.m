function pv = pv_launch(s)
    % PV_LAUNCH()
    %   This is a 'launcher' file for PoreView. It is designed to start an
    %   instance of PoreView in a specified folder and to give it the 
    %   keyboard callback behavior you want.
    
    global fcurrentSigNoSpikes
    
    % this sets the default directory for File->Open
    if nargin < 1
        s = '/Users/Stephen/Documents/Stephen/Research/Data/Biopore/';
    end
	pv = PoreView(s);
    
    % variable to hold the ranges we are trimming
    ranges = [];
    
    function keyFn(e)
        
        % do this even if we don't have data loaded
        if strcmp(e.Character,'o')
            % make a nice plotted version of the PoreView window            
            plot_signals(pv);
        end
        
        % do nothing if we don't have data loaded yet
        if isempty(pv.data)
            return
        end
        
        % to figure out what the keys are called, uncomment this line
        %disp(e);
        
        if strcmp(e.Character,'f')
            % ask user what filters they want to add

            str = inputdlg('Enter desired filter and frequency/param:      ','PoreView',1,{'lpb 1000'});
            
            strs = strsplit(str{1});
            
            if numel(strs) < 2
                return
            elseif numel(strs)==3 % band pass
                params(1) = str2double(strs{2});
                params(2) = str2double(strs{3});
                param = 1;
            else
                param = str2double(strs{2});
                if isnan(param) || param <= 0
                    return
                end
            end
            
            switch strs{1}
                case 'lp'
                    filtname = sprintf('Low-pass (%d Hz)', param);
                    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_lp(d,4,param),filtname), pv.data, 'UniformOutput', false);
                    fsigs = fsigs{1};
                case 'lpb'
                    filtname = sprintf('Low-pass Bessel (%d Hz)', param);
                    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_lpb(d,4,param),filtname), pv.data, 'UniformOutput', false);
                    fsigs = fsigs{1};
                case 'hp'
                    filtname = sprintf('High-pass (%d Hz)', param);
                    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_hp(d,4,param),filtname), pv.data, 'UniformOutput', false);
                    fsigs = fsigs{1};
                case 'med'
                    filtname = sprintf('Median (%d pts)', param);
                    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_med(d,param),filtname), pv.data, 'UniformOutput', false);
                    fsigs = fsigs{1};
                case 'band'
                    filtname = sprintf('Band-pass (%d Hz)', params);
                    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_band(d,4,params),filtname), pv.data, 'UniformOutput', false);
                    fsigs = fsigs{1};
                otherwise
                    return
            end            
            
            % and replace original signals with new ones
            % uh trust me on this one
            for i=1:numel(pv.psigs)
                s = pv.psigs(i).sigs;
                s(s<=pv.data(1).nsigs+1) = s(s<=pv.data(1).nsigs+1) + fsigs(1) - 2;
                pv.psigs(i).sigs = s;
            end
            
            pv.refresh();

        elseif strcmp(e.Character,'n')
            % display a noise plot, a la ClampFit

            % if cursors, do those
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = pv.getView();
            end
            % then make a noise plot
            plot_noise(pv.data,tr);

        elseif strcmp(e.Character,'s')
            % select channels (in a fast5 file)

            if ~strcmp(pv.data.ext,'.fast5')
                return
            end

            % pop up input box to enter desired channels
            val = inputdlg('Enter channels to view, separated by commas:      ','PoreView');
            if isempty(val) || isempty(val{1})
                return
            end
            chans = str2num(val{1});
            chans = chans(and(chans >= pv.data.header.minChan, chans <= pv.data.header.maxChan));
            if isempty(chans)
                errordlg('Invalid channel numbers entered!','PoreView');
                return
            end
            % reload with appropriate channels

            pv.data = SignalData(pv.data.filename,'Channels',chans);
            % remove any signals that exist no more
            for i=1:numel(pv.psigs)
                pv.psigs(i).sigs = pv.psigs(i).sigs(pv.psigs(i).sigs <= numel(chans));
                if isempty(pv.psigs(i).sigs)
                    % and set to some default, so no empty panels
                    pv.psigs(i).sigs = 2;
                end
            end
            pv.refresh();

        % Extras - Stephen Fleming ----------------------------------------
        
        elseif strcmp(e.Character,'k')
            % remove a range of points between cursors
            xx = pv.getCursors();
            if isempty(xx)
                % cursors are invisible
                return
            end
            
            % get average of endpoints, in a narrow range around them
            y0s = mean(pv.data.getByTime(xx(1),xx(1)+0.001));
            y1s = mean(pv.data.getByTime(xx(2),xx(2)-0.001));
            % and their average
            yave = mean([y0s; y1s]);
            % and add it to ranges
            ranges(end+1,:) = [xx yave(2:pv.data.nsigs+1)];
            % update virtual signal
            range = pv.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
            % and refresh visible points
            pv.setSignalPanel(1, range(1));
            pv.setSignalPanel(2, range(2));
            pv.refresh();
            % and display some stuff
            fprintf('Removed %f to %f\n',xx(1),xx(2));
            display(ranges)
            
        elseif strcmp(e.Character,'p')
            % plot in a print-worthy way
            % if cursors, do those
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = pv.getView();
            end
            util.doPlot(pv,tr);
        
        elseif strcmp(e.Character,'l')
            % find discrete levels in data
            % if cursors, do those
            tr = pv.getCursors();
%             filter = 1000;
%             sample = 10000;
%             pValue = -10;
%             levels = util.doLevelAnalysis(pv.data, tr, filter, sample, pValue);
            inpt = inputdlg({'Expected level duration (ms)','Tolerable false positives per second'},'Input level characteristics',1,{'430','1e-4'});
            sigs = pv.data.getSignalList(); % list signals currently displayed in panel 1
            sig = sigs{pv.psigs(1).sigs(end)}; % get the name of the top signal in panel 1
            nums = regexp(sig, '\d+', 'match');
            filter = str2double(nums{1}); % get the filter frequency from the signal name
            if filter==0
                filter = 10000;
            end
            display(['Detected a filter of ' num2str(filter) 'Hz'])
            current = util.downsample_pointwise(pv.data, pv.psigs(1).sigs(end), tr, min(filter*5*diff(tr),10e6));
            t = linspace(tr(1),tr(2),numel(current));
            data = [t', current'];
            levels = karplus_levels(data, 1/(str2double(inpt{1})/1000), str2double(inpt{2}), filter);
            assignin('base','levels',levels);
            % plot on poreview
            plot(pv.psigs(1).axes, cell2mat(cellfun(@(x) [x.start_time, x.end_time], levels, 'uniformoutput', false))', ...
                 cell2mat(cellfun(@(x) ones(1,2)*x.current_median, levels, 'uniformoutput', false))', 'Color', 'r','LineWidth',2);
            plot(pv.psigs(1).axes, cellfun(@(x) x.start_time, levels)', cellfun(@(x) x.current_median, levels)', 'go');
            assignin('base','levels',levels); % assign variable to workspace
            
        elseif strcmp(e.Character,'m')
            % display the mean
            % between cursors
            tr = pv.getCursors();
%             current = util.downsample_pointwise(pv.data,4,tr,1000)*1000;
%             current = pv.data.get(round(tr(1)/pv.data.si):round(tr(2)/pv.data.si), ...
%                 pv.psigs(1).sigs)*1000;
            current = util.downsample_pointwise(pv.data,pv.psigs(1).sigs,tr,1e6)*1000;
            V = mean(util.downsample_pointwise(pv.data,3,tr,1000));
            m = mean(current);
            stdev = std(current);
            display(' ')
            display(['Mean =  ' num2str(m,5) ' ' char(177) ' ' num2str(stdev,3) ' pA'])
            display(['Conductance =  ' num2str(m/V,3) ' nS'])
            display(['Voltage =  ' num2str(V,5) ' mV'])
            
        elseif strcmp(e.Character,'x')
            % display conductance
            display('Showing conductance')
            % create a conductance channel for all signals that are current
            filtname = 'Conductance (nS)';
            fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_conductance(d,3), ...
                filtname, [2 3]), pv.data, 'UniformOutput', false);
            if pv.data.nvsigs>0 % have to do it this way for order of operations not to screw things up
                % go through each virtual signal and make a new one for
                % conductance, that uses only the original data, and that
                % names it appropriately
                for i = 1:(pv.data.nvsigs/2)
                    % if this is conductance already, skip it
                    if strfind(pv.data.vnames{i},'Conductance')
                        continue;
                    end
                    fun = pv.data.vfuns{i};
                    arrayfun(@(x) x.addVirtualSignal(@(d) filt_conductance(fun(d),3), ...
                        sprintf('%s (%s)',filtname,pv.data.vnames{i}), [2 3]), pv.data, 'UniformOutput', false);
                end
            end
            fcondSig = fsigs{1};
            pv.setSignalPanel(1, fcondSig(1)); % show it
            pv.autoscaleY();
        
        elseif strcmp(e.Character,'h')
            % do a histogram
            % between cursors
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = pv.getView();
            end
            util.doHistogram(pv.data,tr);
            
        elseif strcmp(e.Character,'q')
            % play audio
            tr = pv.getCursors();
            r = round(1/(20000*pv.data.si));
            display(['Playing back audio: [' num2str(tr) '], sampling rate ' num2str(1/pv.data.si/r)]);
            sig = listdlg('PromptString','Select a signal:',...
                'SelectionMode','single',...
                'ListString',pv.data.getSignalList);
            n = 1e6;
            red = pv.data.getViewData(tr);
            m = max(red(:,sig));
            for i=round(tr(1)/pv.data.si):n:round(tr(2)/pv.data.si)
                %display(['Playing: ' num2str(round((i*pv.data.si))) ' to ' num2str(round(min(i+n,round(tr(2)/pv.data.si))*pv.data.si))]);
                d = pv.data.get(i:min(i+n,round(tr(2)/pv.data.si)),sig);
                d = d/m; % normalize to 1
                d = decimate(d,r);
                sound = audioplayer(d,1/pv.data.si/r);
                playblocking(sound);
            end
            
        elseif strcmp(e.Character,'r')
            % perform a low-voltage diffusion analysis
            tr = pv.getCursors();
            util.doLowVoltageDiffusionAnalysis(pv,tr);
            
        elseif strcmp(e.Character,'g')
            % generate 'spike' information from an ideal capacitive current spike
            % if cursors, do those
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = [];
                display('No region selected!')
            else
                display(['Generating spike template from [' num2str(tr) ']'])
                out = remove_spikes(pv.data,2,tr,[],0,1);
                spike = out(2:end);
                assignin('base','spike',spike);
                spikeIndex = out(1);
                assignin('base','spikeIndex',spikeIndex);
                file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                    pv.data.filename(65:68) pv.data.filename(70:71) pv.data.filename(73:74) '/' ...
                    pv.data.filename(76:end-4) '_spike.mat'];
                save(file,'spike','spikeIndex');
                display('Saved spike template.')
                
                % remove capacitive current spikes throughout data, creating new virtual signal
                %display('Removing capacitive current spikes.')
                %fcurrentSigNoSpikes = pv.data.addVirtualSignal(@(d) repmat(remove_spikes([d(:,1),d(:,2),d(:,3)],3,[],spike,spikeIndex,0),[1 3]), 'Spikes removed', [1 3 2]);
            end
            
        elseif strcmp(e.Character,'z')
            % have cursors select all
            tr = [pv.data.tstart, pv.data.tend];
            pv.setCursors(tr);
            pv.setView(tr);
            
        elseif strcmp(e.Character,'t')
            % perform analysis of blockage data gathered at a series of
            % voltages
            tr = pv.getCursors();
            util.doBlockageAnalysis(pv.data,tr);
            
        elseif strcmp(e.Character,'e')
            % put cursors at edges of the event that the cursors are
            % currently inside
            tr = pv.getCursors();
            util.doFindEventEdges(pv,tr);
            
        elseif strcmp(e.Character,'v')
            tr = pv.getCursors();
            filter = 100;
            current = util.downsample_pointwise(pv.data, pv.psigs(1).sigs(end), tr, min(filter*5*diff(tr),10e6));
            t = linspace(tr(1),tr(2),numel(current));
            data = [t', current'];
            clear t current
            levels = karplus_levels(data, 1/(1600/1000), 1e-10, filter);
            clear data
            [k, pulses] = pulse_correlation(pv.data,tr,4,levels);
            disp(['kappa = ' num2str(k)])
            
        end
        
    end
    
    % and set our all-important keyboard callback
    pv.setKeyboardCallback(@keyFn);
end



