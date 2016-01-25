classdef util
    
    methods (Static = true)
        
        function doPlot(pv,tr)
            
            t = inputdlg('Enter title:','Input');
            display(['Plotting time interval [' num2str(tr) ']'])
            plot_pretty(pv.data,tr,1000,[5,3],t);
            %plot_pretty(pv.data,tr,'none',[4,3],t);
            %plot_pretty(pv.data,tr,2000,[2,3,6],t);
            
        end
        
        function doFindEvents(pv,minLength)
            
            % find all events (gaps in open pore) that are longer than
            % minLength
            % assumes open pore is recorded in the file
            
            % first find the open pore current by doing a histogram
            tr = [pv.data.tstart, pv.data.tend];
            del = 0.1; % in pA
            raw = pv.data.getViewData(tr);
            lowCurrent = min(raw(:,2)*1000);
            highCurrent = max(raw(:,2)*1000);
            xcurrent = lowCurrent-5:del:highCurrent+5;
            hcurrent = hist(raw(:,2)*1000,xcurrent);
            i1 = find(hcurrent>50,1,'last');
            i2 = find(hcurrent(1:i1)<50,1,'last');
            g = figure(2);
            set(g,'position',[-987         136         908        1067])
            subplot(4,1,1)
            bar(xcurrent,hcurrent);
            hold on
            bar(xcurrent(i2:i1),hcurrent(i2:i1),'r','EdgeColor','r')
            xlim([-20 xcurrent(i1)+20]);
            xlabel('Current (pA)')
            ylabel('Data points')
            set(gca,'fontsize',16)
            % get mean open pore current
            mean_open = sum(hcurrent(i2:i1)*del.*xcurrent(i2:i1))/sum(hcurrent(i2:i1)*del);
            plot(mean_open,3*max(hcurrent(i2:i1)),'r*');
            
        end
        
        function doFindEventEdges(pv,tr)
            
            raw1 = abs(pv.data.getViewData([max(0,tr(1)-500),tr(1)]));
            raw2 = abs(pv.data.getViewData([tr(2),min(pv.data.tend,tr(2)+500)]));
            threshold = 0.4;
            rawInd1 = find(raw1(:,2)>threshold,1,'last');
            if isempty(rawInd1)
                trange(1) = 0;
            else
                trange(1) = pv.data.findPrev(@(x) abs(x(:,2))>threshold, (raw1(rawInd1,1)+0.05)/pv.data.si) * pv.data.si;
                if trange(1)<0
                    trange(1) = raw1(rawInd1,1); % default to the raw guess if we get something wrong
                end
            end
            rawInd2 = find(raw2(:,2)>threshold,1,'first');
            if isempty(rawInd2)
                trange(2) = pv.data.tend;
            else
                trange(2) = pv.data.findNext(@(x) abs(x(:,2))>threshold, (raw2(rawInd2,1)-0.05)/pv.data.si) * pv.data.si;
            end
            pv.setCursors(trange);
            display(['Time = ' num2str(trange(2)-trange(1),5)])
            
        end
        
        function doBlockageAnalysis(pv,tr)
            
            display(['Voltage-dependent current blockage anaylsis of data [' num2str(tr) ']'])
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                pv.data.filename(65:68) pv.data.filename(70:71) pv.data.filename(73:74) '/' pv.data.filename(76:end-4) ...
                '_blockage_info.mat'];
            cr = [0, 1, 2]; % sepcify conductance range of one blockage [low, high, open]
            [blocks, open_pore] = blockage_analysis(pv.data,tr,cr);
            assignin('base','blocks',blocks);
            assignin('base','open_pore',open_pore);
            answer = input('Save data? (y/n): ','s');
            if strcmp(answer,'y')
                save(file,'blocks','open_pore');
                fprintf('\n')
                display(['Saved data as ' file])
            end
            
        end
        
        function doLowVoltageDiffusionAnalysis(pv,tr)
            
            display(['Low-voltage diffusion anaylsis of data [' num2str(tr) ']'])
            channel = fcurrentSigNoSpikes(2);
            %fcondSig = pv.data.addVirtualSignal(@(d) repmat(d(:,2)*1000./max(abs(d(:,3)),5).*sign(d(:,3)),[1 2]),'Conductance (nS)',[channel,3]); % signal 6
            fcondSig = pv.data.addVirtualSignal(@(d) repmat(d(:,2)*1000./d(:,3),[1 2]),'Conductance (nS)',[channel,3]); % signal 6
            pv.setSignalPanel(1,fcondSig(1));
            pv.psigs(1).setY([-0.5 3]);
            pv.refresh();
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                pv.data.filename(65:68) pv.data.filename(70:71) pv.data.filename(73:74) '/' pv.data.filename(76:end-4) ...
                '_event_info.mat'];
            captureVoltage = 160; % in mV
            holdingVoltageRange = [25 75]; % in mV
            conductanceCutoff = 1.5; % in nS, if captured molecule conductance is greater than this, it's thrown out
            events = low_voltage_diffusion(pv.data,tr,captureVoltage,holdingVoltageRange,conductanceCutoff,fcondSig);
            assignin('base','events',events);
            answer = input('Save data? (y/n): ','s');
            if strcmp(answer,'y')
                save(file,'events');
                fprintf('\n')
                display(['Saved data as ' file])
            end
            
        end
        
        function doCurrentHistogram(pv,tr)
            
            display(['Creating histogram from [' num2str(tr) ']'])
            if ~exist('fcurrentSig','var')
                filter = 10000; % Hz
                rangeEdited = pv.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
                fcurrentSig = pv.data.addVirtualSignal(@(d) filt_lp(d,4,filter)*1e3,'Low-pass',rangeEdited(1)); % signal 5
                fcondSig = pv.data.addVirtualSignal(@(d) repmat(d(:,2)./d(:,3),[1 2]),'Conductance (nS)',[fcurrentSig,rangeEdited(2)]); % signal 6
            end
            channels = [fcurrentSig,fcondSig];
            histogram_pv(pv.data,tr,filter,channels);
            
        end
        
        function doLevelAnalysis(pv,tr,finalFrequency,pValue)
            
            display('Finding discrete levels:')
            discreteData = find_discrete_levels(pv.data, 4, tr, finalFrequency, pValue);

            h = figure(2);
            clf(2)
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
            Vdata = pv.data.getViewData([tr(1),(tr(1)+0.001)]);
            V = round(mean(Vdata(:,3))/5)*5;
            assignin('base','V',V);
            name = [pv.data.filename(65:68) '\_' pv.data.filename(70:71) '\_' pv.data.filename(73:74) '\_' pv.data.filename(76:end-4)];
            assignin('base','name',name);
            assignin('base','tr',tr);
            figure(2)
            h = get(gca,'Title');
            title = get(h,'String');
            plot_squiggles(discreteData, name, tr, title); % plot the level information
            %plot_level_duration(discreteData, name, tr, title); % plot step duration distribution
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
            
        end
        
    end
    
end