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

            str = inputdlg('Enter desired filter and frequency/param:      ','PoreView',1,{'lp 10000'});
            
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
                    fsigs = pv.data.addVirtualSignal(@(d) filt_lp(d,4,param),filtname);
                case 'lpb'
                    filtname = sprintf('Low-pass Bessel (%d Hz)', param);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_lp(d,4,param),filtname);
                case 'hp'
                    filtname = sprintf('High-pass (%d Hz)', param);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_hp(d,4,param),filtname);
                case 'med'
                    filtname = sprintf('Median (%d pts)', param);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_med(d,param),filtname);
                case 'band'
                    filtname = sprintf('Band-pass (%d Hz)', params);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_band(d,4,params),filtname);
                otherwise
                    return
            end            
            
            % and replace original signals with new ones
            % uh trust me on this one
            for i=1:numel(pv.psigs)
                s = pv.psigs(i).sigs;
                s(s<=pv.data.nsigs+1) = s(s<=pv.data.nsigs+1) + fsigs(1) - 2;
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
            xlim = pv.getCursors();
            if isempty(xlim)
                % cursors are invisible
                return
            end
            
            % get average of endpoints, in a narrow range around them
            y0s = mean(pv.data.getByTime(xlim(1),xlim(1)+0.001));
            y1s = mean(pv.data.getByTime(xlim(2),xlim(2)-0.001));
            % and their average
            yave = mean([y0s; y1s]);
            % and add it to ranges
            ranges(end+1,:) = [xlim yave(2:pv.data.nsigs+1)];
            % update virtual signal
            range = pv.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
            % and refresh visible points
            pv.setSignalPanel(1, range(1));
            pv.setSignalPanel(2, range(2));
            pv.refresh();
            % and display some stuff
            fprintf('Removed %f to %f\n',xlim(1),xlim(2));
            display(ranges)
            
        elseif strcmp(e.Character,'p')
            % plot in a print-worthy way
            % if cursors, do those
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = pv.getView();
            end
            t = inputdlg('Enter title:','Input');
            display(['Plotting time interval [' num2str(tr) ']'])
            %plot_pretty(pv.data,tr,1000,[2,3],t);
            plot_pretty(pv.data,tr,1000,2,t);
            %plot_pretty(pv.data,tr,2000,[2,3,6],t);
        
        elseif strcmp(e.Character,'l')
            % find discrete levels in data
            % if cursors, do those
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = pv.getView();
            end
            t = inputdlg('Enter title:','Input');
            % add filtered signals
            filterFreq = 100;
            fs = 1000;
            if ~exist('f2','var')
                f1 = pv.data.addVirtualSignal(@(d) filt_lpb(d,4,filterFreq),'Low-pass Bessel 100Hz',[1,2]);
                f2 = pv.data.addVirtualSignal(@(d) filt_med(d,51),'Median filtered',f1);
                assignin('base','f2',f2);
            end
            display(['Data filtered at ' num2str(filterFreq) 'Hz and downsampled to ' num2str(fs) 'Hz.'])
            [discreteData, V] = analyze_level_data2(pv.data,tr,f2,fs,t);
            assignin('base','discreteData',discreteData); % assign variable to workspace
            assignin('base','V',V);
            name = [pv.data.filename(65:68) '\_' pv.data.filename(70:71) '\_' pv.data.filename(73:74) '\_' pv.data.filename(76:end-4)];
            assignin('base','name',name);
            assignin('base','tr',tr);
            figure(2)
            h = get(gca,'Title');
            title = get(h,'String');
            plot_squiggles(discreteData, name, tr, title); % plot the level information
            plot_level_duration(discreteData, name, tr, title); % plot step duration distribution
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                pv.data.filename(65:68) pv.data.filename(70:71) pv.data.filename(73:74) '/' pv.data.filename(76:end-4) ...
                '_discreteData_' num2str(round(tr(1))) '.mat'];
            
            % Add pulsing data if we have it
            if pv.data.nsigs > 2 % we have pulsing channel
                [pulses, candidates, distances] = pulse_analysis(pv.data, tr, discreteData, [0.005, -0.002]);
                figure(2)
                line(repmat(candidates',1,2)',repmat(get(gca,'ylim'),numel(candidates),1)','Color',[0 1 0],'LineStyle','-') % vertical lines
                annotation('textbox', [0.8 0.85 0 0], 'String', [name ' ' num2str(numel(candidates)) '/' num2str(numel(pulses))], 'FontSize', 20);
                pulse.times = pulses;
                pulse.candidates = candidates;
                pulse.distances = distances;
                assignin('base','pulse',pulse);
                answer = input('Save data? (y/n): ','s');
                if strcmp(answer,'y')
                    save(file,'discreteData','V','name','tr','pulse')
                    display(['Saved data as ' file])
                end
            else
                answer = input('Save data? (y/n): ','s');
                if strcmp(answer,'y')
                    save(file,'discreteData','V','name','tr')
                    display(['Saved data as ' file])
                end
            end
            
        elseif strcmp(e.Character,'m')
            % display the mean
            % between cursors
            tr = pv.getCursors();
            [data,~] = downsample(pv.data,2,tr,1000,10000);
            V = mean(downsample(pv.data,3,tr,100,1000));
            m = mean(data,2);
            s = std(data,1,2);
            display(['Mean =  ' num2str(m,3) ' ' char(177) ' ' num2str(s,3) ' pA'])
            display(['Conductance =  ' num2str(m/V,3) ' nS'])
            display(['Voltage =  ' num2str(V) ' mV'])
            
        elseif strcmp(e.Character,'x')
            % display conductance
            display('Showing conductance')
            %if ~exist('fcurrentSig','var')
                filter = 1000; % Hz
                fcurrentSig = pv.data.addVirtualSignal(@(d) filt_lp(d,4,filter)*1e3,'Low-pass',2); % signal 5
                fcondSig = pv.data.addVirtualSignal(@(d) repmat(d(:,2)./d(:,3),[1 2]),'Conductance (nS)',[fcurrentSig,3]); % signal 6
            %end
            pv.setSignalPanel(3, fcondSig(2)); % show it
            pv.autoscaleY();
        
        elseif strcmp(e.Character,'h')
            % do a current histogram
            % between cursors
            tr = pv.getCursors();
            display(['Creating histogram from [' num2str(tr) ']'])
            if ~exist('fcurrentSig','var')
                filter = 10000; % Hz
                rangeEdited = pv.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
                fcurrentSig = pv.data.addVirtualSignal(@(d) filt_lp(d,4,filter)*1e3,'Low-pass',rangeEdited(1)); % signal 5
                fcondSig = pv.data.addVirtualSignal(@(d) repmat(d(:,2)./d(:,3),[1 2]),'Conductance (nS)',[fcurrentSig,rangeEdited(2)]); % signal 6
            end
            channels = [fcurrentSig,fcondSig];
            histogram_pv(pv.data,tr,filter,channels);
            
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
            display(['Low-voltage diffusion anaylsis of data [' num2str(tr) ']'])
            channel = fcurrentSigNoSpikes(2);
            fcondSig = pv.data.addVirtualSignal(@(d) repmat(d(:,2)*1000./d(:,3),[1 2]),'Conductance (nS)',[channel,3]); % signal 6
            pv.setSignalPanel(1,fcondSig(1));
            pv.refresh();
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                pv.data.filename(65:68) pv.data.filename(70:71) pv.data.filename(73:74) '/' pv.data.filename(76:end-4) ...
                '_event_info.mat'];
            captureVoltage = 160; % in mV
            holdingVoltageRange = [30 75]; % in mV
            conductanceCutoff = 4; % in nS, if captured molecule conductance is greater than this, it's thrown out
            events = low_voltage_diffusion(pv.data,tr,captureVoltage,holdingVoltageRange,conductanceCutoff,fcondSig);
            assignin('base','events',events);
            answer = input('Save data? (y/n): ','s');
            if strcmp(answer,'y')
                save(file,'events');
                fprintf('\n')
                display(['Saved data as ' file])
            end
            
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
                
                % remove capacitive current spikes throughout data, creating new virtual signal
                display('Removing capacitive current spikes.')
                fcurrentSigNoSpikes = pv.data.addVirtualSignal(@(d) repmat(remove_spikes([d(:,1),d(:,2),d(:,3)],3,[],spike,spikeIndex,0),[1 3]), 'Spikes removed', [1 3 2]);
            end
            
        elseif strcmp(e.Character,'z')
            % have cursors select all
            tr = [pv.data.tstart, pv.data.tend];
            pv.setCursors(tr);
            pv.setView(tr);
            
        elseif strcmp(e.Character,'w')
            % automatic selection of level transitions
            % filters
            pCutoff = input('Type the desired p-value cutoff (e.g. 1e-80): ');
            display('Applying filters...')
            filter = 5; % in Hz
            lp = pv.data.addVirtualSignal(@(d) filt_lpb(d,8,filter), ['Low pass ' num2str(filter) 'Hz'], 2);
            med = pv.data.addVirtualSignal(@(d) filt_med(d,71), 'Median 71pt', lp);
            % find the levels
            levels = pv.data.addVirtualSignal(@(d) filt_discrete_levels(d, pCutoff), 'Levels', med);
            pv.setSignalPanel(1,[med, levels]);
            pv.refresh();
            display('Done.')
            
        elseif strcmp(e.Character,'t')
            % perform analysis of blockage data gathered at a series of
            % voltages
            tr = pv.getCursors();
            display(['Voltage-dependent current blockage anaylsis of data [' num2str(tr) ']'])
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                pv.data.filename(65:68) pv.data.filename(70:71) pv.data.filename(73:74) '/' pv.data.filename(76:end-4) ...
                '_blockage_info.mat'];
            blocks = blockage_analysis(pv.data,tr);
            assignin('base','blocks',blocks);
            answer = input('Save data? (y/n): ','s');
            if strcmp(answer,'y')
                save(file,'blocks');
                fprintf('\n')
                display(['Saved data as ' file])
            end
            
        end
            
    end

    % and set our all-important keyboard callback
    pv.setKeyboardCallback(@keyFn);
end

