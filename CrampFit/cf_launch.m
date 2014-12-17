function cf = cf_launch()
    % CF_LAUNCH()
    %   This is a 'launcher' file for CrampFit. It is designed to start an
    %   instance of CrampFit in a specified folder and to give it the 
    %   keyboard callback behavior you want.

    %% this sets the default directory for File->Open
	cf = CrampFit('/Users/Stephen/Documents/Stephen/Research/Data/Biopore/20141021/');
    
    % variable to hold the ranges we are trimming
    ranges = [];
    
    function keyFn(e)
        % do nothing if we don't have data loaded yet
        if isempty(cf.data)
            return
        end
        
        % to figure out what the keys are called, uncomment this line
        %disp(e);
        
        if strcmp(e.Character,'k')
            % remove a range of points between cursors
            xlim = cf.getCursors();
            if isempty(xlim)
                % cursors are invisible
                return
            end
            
            % get average of endpoints, in a narrow range around them
            y0s = mean(cf.data.getByTime(xlim(1),xlim(1)+0.001));
            y1s = mean(cf.data.getByTime(xlim(2),xlim(2)-0.001));
            % and their average
            yave = mean([y0s; y1s]);
            % and add it to ranges
            ranges(end+1,:) = [xlim yave(2:cf.data.nsigs+1)];
            % update virtual signal
            range = cf.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
            % and refresh visible points
            cf.setSignalPanel(1, range(1));
            cf.setSignalPanel(2, range(2));
            cf.refresh();
            % and display some stuff
            fprintf('Removed %f to %f\n',xlim(1),xlim(2));
            display(ranges)
            
        elseif strcmp(e.Character,'f')
            % create the requisite virtual signals
            
            % subselected data filter
            f_rm = cf.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
            % high pass acts on subselected data
            %f_hp = cf.data.addVirtualSignal(@(d) filt_hp(d,4,200),'High-pass',f_rm);
            % tell median to act on high-passed data
            %f_med = cf.data.addVirtualSignal(@(d) filt_med(d,15),'Median',f_hp);
            f_med = cf.data.addVirtualSignal(@(d) filt_med(d,50),'Median',f_rm);
            
            % also set which signals to draw in each panel, you can play
            % with this all you like
            %cf.setSignalPanel(1, f_rm(1));
            cf.setSignalPanel(1, f_med);
            
            % draw both median-filtered panels
            %cf.addSignalPanel(f_med);

            disp('Filters added')
        
        elseif strcmp(e.Character,'n')
            % display a noise plot!
            
            % if cursors, do those
            tr = cf.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = cf.getView();
            end
            % then make a noise plot
            plot_noise(cf.data,tr);
        
        elseif strcmp(e.Character,'p')
            % plot in a print-worthy way
            % if cursors, do those
            tr = cf.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = cf.getView();
            end
            display(['Plotting time interval [' num2str(tr) ']'])
            %plot_pretty(cf.data,tr,5000,[2,3]);
            plot_pretty(cf.data,tr,5000,2);
        
        elseif strcmp(e.Character,'l')
            % find discrete levels in data
            % if cursors, do those
            tr = cf.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = cf.getView();
            end
            [discreteData, V] = analyze_level_data(cf.data,tr);
            assignin('base','discreteData',discreteData); % assign variable to workspace
            assignin('base','V',V);
            name = [cf.data.filename(65:68) '\_' cf.data.filename(70:71) '\_' cf.data.filename(73:74) '\_' cf.data.filename(76:end-4)];
            assignin('base','name',name);
            assignin('base','tr',tr);
            plot_squiggles(discreteData, V, name, tr); % plot the level information
            plot_level_duration(discreteData, V, name, tr); % plot step duration distribution
            file = ['/Users/Stephen/Documents/Stephen/Research/Analysis/Biopore/' ...
                cf.data.filename(65:68) cf.data.filename(70:71) cf.data.filename(73:74) '/' cf.data.filename(76:end-4) ...
                '_discreteData_' num2str(round(tr(1))) '.mat'];
            
            % Add pulsing data if we have it
            if cf.data.nsigs > 2 % we have pulsing channel
                [pulses, candidates, distances] = pulse_analysis(cf.data, tr, discreteData, [0.005, -0.002]);
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
            tr = cf.getCursors();
            [data,~] = downsample(cf.data,3,tr,1000,10000);
            m = mean(data,2);
            s = std(data,1,2);
            display(['Mean =  ' num2str(m,3) ' ' char(177) ' ' num2str(s,3)])
            
        elseif strcmp(e.Character,'s')
            % display conductance
            display('Showing conductance')
            if ~exist('fcurrentSig','var')
                filter = 1000; % Hz
                fcurrentSig = cf.data.addVirtualSignal(@(d) filt_lp(d,4,filter)*1e3,'Low-pass',2); % signal 5
                fcondSig = cf.data.addVirtualSignal(@(d) repmat(d(:,2)./d(:,3),[1 2]),'Conductance (nS)',[fcurrentSig,3]); % signal 6
            end
            cf.setSignalPanel(3, fcondSig(1)); % show it
            cf.autoscaleY();
        
        elseif strcmp(e.Character,'h')
            % do a current histogram
            % between cursors
            tr = cf.getCursors();
            display(['Creating histogram from [' num2str(tr) ']'])
            filter = 1000; % Hz
            if ~exist('fcurrentSig','var')
                filter = 1000; % Hz
                rangeEdited = cf.data.addVirtualSignal(@(d) filt_rmrange(d,ranges),'Range-edited');
                fcurrentSig = cf.data.addVirtualSignal(@(d) filt_lp(d,4,filter)*1e3,'Low-pass',rangeEdited(1)); % signal 5
                fcondSig = cf.data.addVirtualSignal(@(d) repmat(d(:,2)./d(:,3),[1 2]),'Conductance (nS)',[fcurrentSig,rangeEdited(2)]); % signal 6
            end
            channels = [fcurrentSig,fcondSig];
            histogram(cf.data,tr,filter,channels);
            
        end
        
    end

    % and set our all-important keyboard callback
    cf.setKeyboardCallback(@keyFn);
end

