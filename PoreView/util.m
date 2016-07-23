classdef util
    
    methods (Static = true)
        
        function doPlot(pv,tr)
            
            t = inputdlg('Enter title:','Input');
            display(['Plotting time interval [' num2str(tr) ']'])
            plot_pretty(pv.data,tr,1000,[5,3],t);
            %plot_pretty(sigdata,tr,'none',[4,3],t);
            %plot_pretty(sigdata,tr,2000,[2,3,6],t);
            
        end
        
        function events = doFindEvents(sigdata, minLength)
            
            % find all events (gaps in open pore) that are longer than
            % minLength
            % assumes open pore is recorded in the file
            
            % first find the open pore current by doing a histogram
            tr = [sigdata.tstart, sigdata.tend];
            del = 0.01; % in nS
            raw = sigdata.getViewData(tr);
            voltage = medfilt1(raw(:,3),10); % get rid of weird spikes (artifact of SignalData min/max downsampling)!!
            if (sigdata.nsigs > 2) % we have pulses
                display('downsampling...')
                % this is much faster than doing real median downsampling
                current = util.downsample_pointwise(sigdata, 2, tr, numel(voltage))'*1000; % do this to avoid pulsing errors
                current = medfilt1(current,30);
            else
                current = medfilt1(raw(:,2),30)*1000;
            end
            % let's just work with positive currents and voltages...
            current = abs(current);
            voltage = abs(voltage);
            conductance = current ./ voltage;
            % get mean open pore conductance
            V = 120; % lowest sliding voltage
            dt = raw(2,1)-raw(1,1);
            clear raw;
            lowCond = max(0.1,min(conductance));
            highCond = min(5,max(conductance));
            xcond = lowCond-0.1:del:highCond+0.1;
            hcond = hist(conductance(voltage>5 & conductance<highCond),xcond);
            [~,ind] = max(hcond(hcond>1.9));
            mean_open = xcond(ind);
            if isempty(mean_open) || (mean_open < 1)
                mean_open = 2.1;
            end
            % find the periods with current blocked between open levels,
            % and voltage high
            blockEndInd = 1;
            i = 1;
            events = cell(0);
            while ~isempty(blockEndInd)
                % find a blockage
                goodEvent = false;
                endedManually = false;
                goesBeyondFile = false;
                voltageHighInd = find(voltage(blockEndInd:end) > 0.9*V, 1, 'first')+blockEndInd;
                if (isempty(voltageHighInd))
                    blockEndInd = [];
                    continue;
                end
                blockStartInd = find(conductance(voltageHighInd:end) < 0.8*mean_open, 1, 'first')+voltageHighInd;
                if (isempty(blockStartInd))
                    blockEndInd = [];
                else
                    blockEndInd = find(conductance(blockStartInd+2:end) > 0.9*mean_open, 1, 'first')+blockStartInd+2;
                end
                % check if the voltage goes below 40mV before the end, if so,
                % reject the event or change the endpoint
                v_lower_lim = 40; % only this low because of weird filter effect on ViewData...
                if (~isempty(blockEndInd))
                    voltageDropInd = find(voltage(blockStartInd+10:blockEndInd)<v_lower_lim,1,'first')+blockStartInd+10;
                    if (~isempty(voltageDropInd) && voltageDropInd < blockEndInd) % if the voltage does in fact drop before the end
                        blockEndInd = voltageDropInd;
                        endedManually = true;
                    end
                else
                    voltageDropInd = find(voltage(blockStartInd+10:end)<v_lower_lim,1,'first')+blockStartInd+10;
                    blockEndInd = voltageDropInd;
                    if (~isempty(voltageDropInd))
                        endedManually = true;
                    else % blockEndInd and voltageDropInd are both nonexistent... the molecule goes past end of file
                        goesBeyondFile = true;
                    end
                end
                if (~isempty(blockStartInd))
                    % time cutoff
                    if (isempty(blockEndInd) && goesBeyondFile)
                        blockEndInd = numel(current);
                    end
                    duration = dt*(blockEndInd - blockStartInd);
                    if (duration > minLength) % less than cutoff time, not real
                        % check if this could be continuing from previous
                        % file... could be past abasic
                        if (blockStartInd < 10)
                            % if it's from the start, then just look for
                            % molecule sliding, not abasic
                            intoMoleculeInd = find(conductance(blockStartInd:blockEndInd)/mean_open<0.35,1,'first');
                            if (~isempty(intoMoleculeInd) && intoMoleculeInd < blockEndInd-blockStartInd)
                                goodEvent = true;
                            end
                        end
                        % look for abasic start
                        start_cond = mean(conductance(blockStartInd:round(blockStartInd+0.05/dt)));
                        if (start_cond/mean_open > 0.42 && start_cond/mean_open < 0.65)
                            % yes, abasic start, but is there more?
                            intoMoleculeInd = find(conductance(blockStartInd:blockEndInd)/mean_open<0.35,1,'first')+blockStartInd-1;
                            if (~isempty(intoMoleculeInd) && (blockEndInd-intoMoleculeInd)*dt > 1) % must have at least a second of good stuff
                                goodEvent = true;
                            end
                        end
                    end
                end
                % if it's a good event, save that info
                if (goodEvent)
                    events{i}.start_ind = blockStartInd;
                    events{i}.end_ind = blockEndInd;
                    events{i}.start_time = blockStartInd*dt;
                    events{i}.end_time = blockEndInd*dt;
                    events{i}.ended_manually = endedManually;
                    events{i}.continues_past_end_of_file = goesBeyondFile;
                    event_voltage = voltage(blockStartInd:blockEndInd);
                    events{i}.voltage = mode(round(event_voltage(event_voltage>80)*100)/100);
                    indCondBaseline = blockEndInd + find(conductance(blockEndInd:end)>mean_open*0.9,1,'first');
                    indCondBaselineEnd = blockEndInd + find(conductance(blockEndInd:end)<mean_open*0.9,1,'first')-1;
                    if ~isempty(indCondBaseline)
                        indVoltageSame = indCondBaseline + find((voltage(indCondBaseline:end)>events{i}.voltage*1.05 | voltage(indCondBaseline:end)<events{i}.voltage*0.95),1,'first');
                        events{i}.open_pore_current = mean(current(indCondBaseline:min([indCondBaselineEnd,indCondBaseline+500,indVoltageSame])));
                    else
                        events{i}.open_pore_current = [];
                    end
                    i = i + 1;
                    display('real molecule found.');
                end
            end
            
            % find the exact starts and ends based on the rough indices
            display('refinement.');
            
            % events were located, go through events
            for j = 1:numel(events)
                % start
                % make sure that if voltage changes, we start at the right spot
                initial = events{j}.start_ind;
                events{j}.start_ind = events{j}.start_ind + find(round(voltage(events{j}.start_ind:events{j}.end_ind))==round(events{j}.voltage),1,'first') - 1;
                inds = events{j}.start_ind + [-5, 5];
                inds = min(max([0,0], dt*inds / sigdata.si),[sigdata.ndata, sigdata.ndata]); % conversion from raw indices to data indices
                full_current = abs(sigdata.get(inds,2)*1000); % current in pA, complete data set
                full_voltage = abs(sigdata.get(inds,3)); % voltage in mV
                full_cond = full_current ./ full_voltage;
                events{j}.start_ind = find(full_cond<0.8*mean_open & round(full_voltage)==round(events{j}.voltage),1,'first')+inds(1)-1;
                if isempty(events{j}.start_ind)
                    events{j}.start_ind = initial;
                    display('strange problem with start ind...')
                end
                % end
                if (~events{j}.ended_manually && ~events{j}.continues_past_end_of_file)
                    inds = events{j}.end_ind + [-5, 5];
                    inds = dt*inds / sigdata.si; % conversion from raw indices to data indices
                    full_current = abs(sigdata.get(inds,2)*1000); % current in pA, complete data set
                    full_voltage = abs(sigdata.get(inds,3)); % voltage in mV
                    full_cond = full_current ./ full_voltage;
                    events{j}.end_ind = find(full_cond>0.9*mean_open,1,'first')+inds(1)-1;
                elseif (~events{j}.ended_manually && events{j}.continues_past_end_of_file)
                    events{j}.end_ind = sigdata.ndata;
                elseif events{j}.ended_manually
                    inds = events{j}.end_ind + [-5, 5];
                    inds = dt*inds / sigdata.si; % conversion from raw indices to data indices
                    full_voltage = abs(sigdata.get(inds,3)); % current in pA, complete data set
                    events{j}.end_ind = find(full_voltage<40,1,'first')+inds(1)-1;
                end
                if isempty(events{j}.end_ind)
                    display('couldn''t find end.  this will cause you problems.')
                    pause();
                end
                events{j}.start_time = events{j}.start_ind * sigdata.si;
                events{j}.end_time = events{j}.end_ind * sigdata.si;
            end
        end
        
        function doFindEventEdges(pv,tr)
            
            sigdata = pv.data;
            raw1 = abs(sigdata.getViewData([max(0,tr(1)-500),tr(1)]));
            raw2 = abs(sigdata.getViewData([tr(2),min(sigdata.tend,tr(2)+500)]));
            threshold = 0.33;
            rawInd1 = find(raw1(:,2)>threshold,1,'last');
            if isempty(rawInd1)
                trange(1) = 0;
            else
                trange(1) = sigdata.findPrev(@(x) abs(x(:,2))>threshold, (raw1(rawInd1,1)+0.05)/sigdata.si) * sigdata.si;
                if trange(1)<0
                    trange(1) = raw1(rawInd1,1); % default to the raw guess if we get something wrong
                end
            end
            rawInd2 = find(raw2(:,2)>threshold,1,'first');
            if isempty(rawInd2)
                trange(2) = sigdata.tend;
            else
                trange(2) = sigdata.findNext(@(x) abs(x(:,2))>threshold, (raw2(rawInd2,1)-0.05)/sigdata.si) * sigdata.si;
            end
            pv.setCursors(trange);
            display(['Time = ' num2str(trange(2)-trange(1),5)])
            
        end
        
        function doBlockageAnalysis(sigdata,tr)
            
            display(['Voltage-dependent current blockage anaylsis of data [' num2str(tr) ']'])
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                sigdata.filename(65:68) sigdata.filename(70:71) sigdata.filename(73:74) '/' sigdata.filename(76:end-4) ...
                '_blockage_info.mat'];
            cr = [0, 1, 2]; % sepcify conductance range of one blockage [low, high, open]
            [blocks, open_pore] = blockage_analysis(sigdata,tr,cr);
            assignin('base','blocks',blocks);
            assignin('base','open_pore',open_pore);
            answer = input('Save data? (y/n): ','s');
            if strcmp(answer,'y')
                save(file,'blocks','open_pore');
                fprintf('\n')
                display(['Saved data as ' file])
            end
            
        end
        
        function doLowVoltageDiffusionAnalysis(sigdata,tr)
            
            display(['Low-voltage diffusion anaylsis of data [' num2str(tr) ']'])
            channel = fcurrentSigNoSpikes(2);
            %fcondSig = sigdata.addVirtualSignal(@(d) repmat(d(:,2)*1000./max(abs(d(:,3)),5).*sign(d(:,3)),[1 2]),'Conductance (nS)',[channel,3]); % signal 6
            fcondSig = sigdata.addVirtualSignal(@(d) repmat(d(:,2)*1000./d(:,3),[1 2]),'Conductance (nS)',[channel,3]); % signal 6
            sigdata.setSignalPanel(1,fcondSig(1));
            sigdata.psigs(1).setY([-0.5 3]);
            sigdata.refresh();
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                sigdata.filename(65:68) sigdata.filename(70:71) sigdata.filename(73:74) '/' sigdata.filename(76:end-4) ...
                '_event_info.mat'];
            captureVoltage = 160; % in mV
            holdingVoltageRange = [25 75]; % in mV
            conductanceCutoff = 1.5; % in nS, if captured molecule conductance is greater than this, it's thrown out
            events = low_voltage_diffusion(sigdata,tr,captureVoltage,holdingVoltageRange,conductanceCutoff,fcondSig);
            assignin('base','events',events);
            answer = input('Save data? (y/n): ','s');
            if strcmp(answer,'y')
                save(file,'events');
                fprintf('\n')
                display(['Saved data as ' file])
            end
            
        end
        
        function doCurrentHistogram(sigdata,tr)
            
            display(['Creating histogram from [' num2str(tr) ']'])
            if ~exist('fcurrentSig','var')
                filter = 10000; % Hz
                fcurrentSig = sigdata.addVirtualSignal(@(d) filt_lp(d,4,filter)*1e3,'Low-pass',2); % signal 5
                fcondSig = sigdata.addVirtualSignal(@(d) repmat(d(:,2)./d(:,3),[1 2]),'Conductance (nS)',[fcurrentSig,2]); % signal 6
            end
            channels = [fcurrentSig,fcondSig];
            histogram_pv(sigdata,tr,filter,channels);
            
        end
        
        function doLevelAnalysis(sigdata,tr,finalFrequency,pvalue)
            
            display('Finding discrete levels:')
            discreteData = find_discrete_levels(sigdata, 2, tr, finalFrequency, pvalue);
            
            h = figure(4);
            clf
            plot(discreteData.time,discreteData.current,'Color',[0.8 0.8 0.8])
            hold on
            line(discreteData.level_timing',(ones(2,1)*discreteData.level_medians'),'Color','k','LineWidth',1)
            cmap = colormap('lines');
            for i = 1:numel(discreteData.level_means)
                i1 = find(discreteData.time >= discreteData.level_timing(i,1),1,'first');
                i2 = find(discreteData.time >= discreteData.level_timing(i,2),1,'first');
                plot(discreteData.time(i1:i2),discreteData.current(i1:i2),'.', ...
                    'MarkerSize',5,'Color',cmap(mod(i-1,size(cmap,1))+1,:))
            end
            ylabel('Current (pA)')
            xlabel('Time (s)')
            xlim([discreteData.time(1) discreteData.time(end)])
            set(h,'Position',[-920 236 700 400]) % size the figure
            set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
            set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
            set(gca,'FontSize',24)
            
            assignin('base','discreteData',discreteData); % assign variable to workspace
            Vdata = sigdata.getViewData([tr(1),(tr(1)+0.001)]);
            V = round(mean(Vdata(:,3))/5)*5;
            assignin('base','V',V);
            name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];
            assignin('base','name',name);
            assignin('base','tr',tr);
            figure(4)
            h = get(gca,'Title');
            title = get(h,'String');
            plot_squiggles(discreteData, name, tr, title); % plot the level information
            %plot_level_duration(discreteData, name, tr, title); % plot step duration distribution
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                sigdata.filename(65:68) sigdata.filename(70:71) sigdata.filename(73:74) '/' sigdata.filename(76:end-4) ...
                '_discreteData_' num2str(round(tr(1))) '.mat'];
            
            % Add pulsing data if we have it
            if sigdata.nsigs > 2 % we have pulsing channel
                [pulses, candidates, distances] = pulse_analysis(sigdata, tr, discreteData, [0.005, -0.002]);
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
            
        end
        
        function r = doFindRegions(sigdata, voltage, conductance, minduration)
            % find all the regions in a data file satisfying certain values
            % of voltage and conductance
            data = sigdata.getViewData([0,min(1740,sigdata.tend)]); % load whole file's downsampled data
            isVoltageOK = and(data(:,3)<voltage+1,data(:,3)>voltage-1);
            isConductanceOK = and(data(:,2)*1000./data(:,3)<conductance, data(:,2)*1000./data(:,3)>0);
            logic = and(isVoltageOK,isConductanceOK);
            % eliminate small islands and blips
            smoothing = 10;
            logic = conv(double(logic),ones(smoothing,1)); % ten in a row needed to change status
            logic = (logic>0);
            dlogic = diff([0; logic; 0]);
            startIndex = find(dlogic > 0);
            endIndex = find(dlogic < 0) - 1;
            dt = data(2,1)-data(1,1);
            duration = (endIndex-startIndex+1) * dt * 1000; % in ms
            if sum(duration>minduration)>0
                rr(:,1) = startIndex(duration>minduration);
                rr(:,2) = endIndex(duration>minduration);
                r = nan(size(rr));
                % refine timings based on real signal
                for i = 1:size(rr,1)
                    % beginning is voltage
                    r(i,1) = sigdata.findPrev(@(d) d(:,3)>voltage+1 | d(:,3)<voltage-1, min((rr(i,1)+smoothing+5),rr(i,2))*dt/sigdata.si);
                    %if r(i,1)*sigdata.si/dt>rr(i,2) % we found something later on
                    %    r(i,1) = nan;
                    %else
                        % end is conductance
                        %r(i,2) = sigdata.findNext(@(d) d(:,2)*1000./d(:,3)>conductance+0.1, max((rr(i,2)-smoothing-10),rr(i,1))*dt/sigdata.si);
                        r(i,2) = sigdata.findNext(@(d) d(:,2)*1000./d(:,3)>conductance+0.1, r(i,1));
                    %end
                end
            else
                r = [];
                return;
            end
            r = r(and(~any(isnan(r),2),(r(:,2)-r(:,1))*sigdata.si>=minduration/1000),:);
        end
        
        function f = getMoleculeFilesAndTimes(mol, varargin)
            
            % determine which files are involved, and how much of each
            % [filenum,  start time, end time;
            %  filenum2, start time, end time;
            %  ...]
            startnum = str2num(mol.start_file(end-7:end-4));
            endnum = str2num(mol.end_file(end-7:end-4));
            f(1,:) = [endnum, 0, mol.end_time];
            for i = endnum-1:-1:startnum
                file = mol.start_file;
                file(end-7:end-4) = sprintf('%04d',i);
                sigdata = SignalData(file);
                f(end+1,:) = [i, 0, sigdata.tend];
            end
            f(end,2) = mol.start_time;
            f = flipud(f);
            
            if numel(varargin) == 1
                % get the time data that corresponds to requested range
                moltimes = cumsum(f(:,3)-f(:,2)); % molecule internal time when it is done with each file
                trange = varargin{1};
                if ~isempty(trange)
                    if numel(trange)==2 && (trange(1)<trange(2)) && diff(trange)<=moltimes(end) % valid input trange
                        firstfile = find(trange(1)<=moltimes,1,'first');
                        lastfile = find(trange(2)<=moltimes,1,'first');
                        if isempty(lastfile)
                            lastfile = firstfile;
                        end
                        whichfilesandtimes = f(firstfile:lastfile,:);
                        moltimes = [0; moltimes];
                        whichfilesandtimes(1,2) = whichfilesandtimes(1,2)+trange(1)-moltimes(firstfile); % start time in first file of interest
                        moltimes = [0; cumsum(whichfilesandtimes(:,3)-whichfilesandtimes(:,2))];
                        whichfilesandtimes(end,3) = whichfilesandtimes(1,2)*(firstfile==lastfile)+diff(trange)-moltimes(end-1); % end is start plus full mol time minus amt of mol before that file
                    else
                        display('Invalid time range.  Displaying full molecule.')
                        trange = [0 moltimes(end)];
                        whichfilesandtimes = f;
                    end
                else
                    whichfilesandtimes = f;
                end
                f = whichfilesandtimes;
            end
            
        end
        
        function d = doLoadMoleculeData(mol, pts, method, varargin)
            % loads current, voltage, and time data for a molecule object
            
            trange = [];
            % check for filter setting
            if numel(varargin)==2
                if numel(varargin{2})==2
                    filter = varargin{1};
                    trange = varargin{2};
                elseif numel(varargin{1})==2
                    filter = varargin{2};
                    trange = varargin{1};
                else
                    filter = varargin{1};
                end
            elseif numel(varargin)==1
                if numel(varargin{1})==2
                    filter = varargin{1};
                elseif numel(varargin{1})==2;
                    trange = varargin{1};
                end
            end
            
            % determine which files are involved, and how much of each
            if isempty(trange)
                filesandtimes = util.getMoleculeFilesAndTimes(mol);
            else
                filesandtimes = util.getMoleculeFilesAndTimes(mol, trange);
            end
            
            % grab data from each file
            t = sum(filesandtimes(:,3)-filesandtimes(:,2)); % total molecule time
            data = [];
            for i = 1:size(filesandtimes,1)
                file = mol.start_file;
                file(end-7:end-4) = sprintf('%04d',filesandtimes(i,1));
                sigdata = SignalData(file);
                
                % filter and smooth
                channel = 2;
                
                % if there are pulses, remove a time range around each to eliminate spikes
                % this is very necessary.  for some reason, spikes thow it off.
                if sigdata.nsigs>2
                    obj.pulses = obj.getPulseTiming();
                    fsigs = sigdata.addVirtualSignal(@(d) filt_rmrange(d,[obj.pulses', obj.pulses'+1.4e-3, nan(size(obj.pulses'))]),'spikes removed');
                    channel = fsigs(1);
                    display('Removing current spikes caused by pulses.')
                end
                
                % downsample
                switch method
                    case 'pointwise'
                        filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter);
                        fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
                        data = cat(1, data, util.downsample_pointwise(sigdata, fsigs(1), filesandtimes(i,2:3), pts*(filesandtimes(i,3)-filesandtimes(i,2))/t)');
                    case 'minmax'
                        filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter);
                        fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
                        data = cat(1, data, util.downsample_minmax(sigdata, fsigs(1), filesandtimes(i,2:3), pts*(filesandtimes(i,3)-filesandtimes(i,2))/t)');
                    case 'median' % takes three times longer than pointwise
                        filtname = sprintf('4-pole Low-pass Bessel (%d Hz)', filter);
                        fsigs = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname);
                        data = cat(1, data, util.downsample_median(sigdata, fsigs(1), filesandtimes(i,2:3), pts*(filesandtimes(i,3)-filesandtimes(i,2))/t)');
                    case 'none' % no downsampling, just grab data
                        d = sigdata.getByTime(filesandtimes(i,2:3));
                        data = cat(1, data, d(:,2));
                end
                
            end
            
            if size(data,2)>size(data,1)
                data = data';
            end
            d = [linspace(0, t, numel(data))', data*1000];
            
        end
        
        function d = doLoadMoleculeViewData(mol, pts, varargin)
            % loads current, voltage, and time data for a molecule object
            
            filesandtimes = util.getMoleculeFilesAndTimes(mol, varargin{:});
            
            % grab data from each file
            t = sum(filesandtimes(:,3)-filesandtimes(:,2)); % total molecule time
            data = [];
            raw = [];
            for i = 1:size(filesandtimes,1)
                file = mol.start_file;
                file(end-7:end-4) = sprintf('%04d',filesandtimes(i,1));
                sigdata = SignalData(file);

                % raw
                raw = cat(1, raw, util.downsample_minmax(sigdata, 'view', filesandtimes(i,2:3), pts*(filesandtimes(i,3)-filesandtimes(i,2))/t)');
            end
            
            d = [linspace(0, t, numel(raw))', raw*1000];
            
        end
        
        function d = doLoadMoleculeViewDataQuick(mol, varargin)
            % loads current, voltage, and time data for a molecule object
            
            filesandtimes = util.getMoleculeFilesAndTimes(mol, varargin{:});
            
            % grab data from each file
            t = sum(filesandtimes(:,3)-filesandtimes(:,2)); % total molecule time
            data = [];
            raw = [];
            for i = 1:size(filesandtimes,1)
                file = mol.start_file;
                file(end-7:end-4) = sprintf('%04d',filesandtimes(i,1));
                sigdata = SignalData(file);

                % raw
                d = sigdata.getViewData([filesandtimes(i,2), filesandtimes(i,3)]);
                raw = cat(1, raw, d(:,2));
            end
            
            d = [linspace(0, t, numel(raw))', raw*1000];
            
        end
        
        function d = downsample_median(sigdata, channel, trange, pts)
            %DOWNSAMPLE_MEDIAN does a median downsampling, returning ABOUT 'pts' points
            start = max(0,ceil(trange(1)/sigdata.si)); % original index
            ending = min(sigdata.ndata,floor(trange(2)/sigdata.si)); % original index
            % downsample data in chunks of 2^22 points
            d = [];
            numpts = 2^22; % data points per chunk
            rep = max(1, round((ending-start) / pts)); % number of original points per downsampled point
            if (rep < 2)
                d = sigdata.get(start:ending,channel);
                return;
            end
            chunks = floor((ending-start)/numpts); % number of full chunks
            fprintf('  0%%\n');
            if chunks ~= 0
                for i = 1:chunks % do chunks of numpts points
                    fulldata = sigdata.get(start+(i-1)*numpts:start+i*numpts-1,channel); % get chunk
                    d = [d, accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@median)'];
                    clear fulldata
                    fprintf('\b\b\b\b\b%3d%%\n',floor(100*i/chunks));
                end
            end
            if mod(pts,numpts)~=0
                fulldata = sigdata.get(start+chunks*numpts:ending,channel); % the last bit that's not a full chunk
                d = [d, accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@median)'];
            end
            fprintf('\b\b\b\b\b');
        end
        
        function d = downsample_pointwise(sigdata, channel, trange, pts)
            %DOWNSAMPLE_POINTWISE does a pointwise downsampling, returning ABOUT 'pts' points
            start = max(0,ceil(trange(1)/sigdata.si)); % original index
            ending = min(sigdata.ndata,floor(trange(2)/sigdata.si)); % original index
            % downsample data in chunks of 2^22
            d = [];
            numpts = 2^22;
            rep = max(1, round((ending-start) / pts)); % number of original points per downsampled point
            if (rep < 2)
                d = sigdata.get(start:ending,channel);
                return;
            end
            chunks = floor((ending-start)/numpts); % number of full chunks
            fprintf('  0%%\n');
            if chunks ~= 0
                for i = 1:chunks % do chunks of numpts points
                    fulldata = sigdata.get(start+(i-1)*numpts:start+i*numpts-1,channel); % get chunk
                    d = [d, downsample(fulldata',rep)];
                    clear fulldata
                    fprintf('\b\b\b\b\b%3d%%\n',floor(100*i/chunks));
                end
            end
            if mod(pts,numpts)~=0
                fulldata = sigdata.get(start+chunks*numpts:ending,channel); % the last bit that's not a full chunk
                d = [d, downsample(fulldata',rep)];
            end
            d = medfilt1(d,10);
            fprintf('\b\b\b\b\b');
        end
        
        function d = downsample_minmax(sigdata, channel, trange, pts)
            %DOWNSAMPLE_MINMAX does a min/max downsampling, returning ABOUT 'pts' points
            
            if strcmp(channel,'view') % only want view data downsampled
                viewdata = sigdata.getViewData(trange);
                viewdata = viewdata(:,2);
                rep = max(1, 2*round(numel(viewdata) / pts)); % number of original points per downsampled point
                d1 = [];
                d2 = [];
                d1 = [d1, accumarray(1+floor((1:numel(viewdata))/rep)',viewdata',[],@max)'];
                d2 = [d2, accumarray(1+floor((1:numel(viewdata))/rep)',viewdata',[],@min)'];
                d = reshape([d1; d2],[1, numel(d1)+numel(d2)]);
                return;
            end
            
            start = max(0,ceil(trange(1)/sigdata.si)); % original index
            ending = min(sigdata.ndata,floor(trange(2)/sigdata.si)); % original index
            % downsample data in chunks of 500000
            d1 = [];
            d2 = [];
            numpts = 2^22;
            rep = max(1, 2*round((ending-start) / pts)); % number of original points per downsampled point
            if (rep < 2)
                d = sigdata.get(start:ending,channel);
                return;
            end
            chunks = floor((ending-start)/numpts); % number of full chunks
            % max
            fprintf('  0%%\n');
            if chunks ~= 0
                for i = 1:chunks % do chunks of numpts points
                    fulldata = sigdata.get(start+(i-1)*numpts:start+i*numpts-1,channel); % get chunk
                    d1 = [d1, accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@max)'];
                    d2 = [d2, accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@min)'];
                    clear fulldata
                    fprintf('\b\b\b\b\b%3d%%\n',floor(100*i/chunks));
                end
            end
            if mod(pts,numpts)~=0
                fulldata = sigdata.get(start+chunks*numpts:ending,channel); % the last bit that's not a full chunk
                d1 = [d1, accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@max)'];
                d2 = [d2, accumarray(1+floor((1:numel(fulldata))/rep)',fulldata',[],@min)'];
            end
            % put them together
            d = reshape([d1; d2],[1, numel(d1)+numel(d2)]);
            fprintf('\b\b\b\b\b');
        end
        
    end
    
end