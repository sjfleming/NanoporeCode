function pv = pv_launch(s)
    % PV_LAUNCH()
    %   This is a 'launcher' file for PoreView. It is designed to start an
    %   instance of PoreView in a specified folder and to give it the 
    %   keyboard callback behavior you want.

    % this sets the default directory for File->Open
    if nargin < 1
        s = 'C:\Minion';
    end
	pv = PoreView(s);
    
    % variable to hold the ranges we are trimming
    ranges = [];
    
    function keyFn(e)
        
        % do this even if we don't have data loaded
        if strcmp(e.Character,'p')
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
            end
            
            param = str2double(strs{2});
            if isnan(param) || param <= 0
                return
            end

            switch strs{1}
                case 'lp'
                    filtname = sprintf('Low-pass (%d Hz)', param);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_lp(d,4,param),filtname);
                case 'hp'
                    filtname = sprintf('High-pass (%d Hz)', param);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_hp(d,4,param),filtname);
                case 'med'
                    filtname = sprintf('Median (%d pts)', param);
                    fsigs = pv.data.addVirtualSignal(@(d) filt_med(d,param),filtname);
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
        end
    end

    % and set our all-important keyboard callback
    pv.setKeyboardCallback(@keyFn);
end

