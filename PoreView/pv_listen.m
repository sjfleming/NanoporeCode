function pv = pv_listen()
    % PV_LISTEN()
    %   This monitors a particular folder/subfolders for new abf/cbf files,
    %   and loads them as it finds them.

    % set the root directory to monitor here
    filedir = 'D:\AxoData\';
	pv = PoreView(filedir);
    
    function keyFn(e)
        % do nothing if we don't have data loaded yet
        if isempty(pv.data)
            return
        end
        
        % to figure out what the keys are called, uncomment this line
        %disp(e);
        
        if strcmp(e.Character,'n')
            % display a noise plot!
            
            % if cursors, do those
            tr = pv.getCursors();
            if isempty(tr)
                % otherwise, do the full view
                tr = pv.getView();
            end
            % then make a noise plot
            plot_noise(pv.data,tr);
        end
    end

    % and set our all-important keyboard callback
    pv.setKeyboardCallback(@keyFn);
    
    % now create a timer to check for new files
    [~,lastfiletime] = getNewestFile(filedir);
    
    function timerFcn(o,~)
        [fn, ftime] = getNewestFile(filedir);
        
        % also stop timer if PoreView closed
        if ~ishandle(pv.fig)
            stop(o);
            delete(o);
            disp('File checking stopped.')
            return;
        end
        
        % do we have a new file that we haven't loaded yet, and also
        % that hasn't been modified in a second or so?
        age = (now - ftime)*24*60*60; % in seconds
        if (ftime > lastfiletime) && (age > 1)
            lastfiletime=ftime;
            pv.loadFile(fn);
        end
    end

    disp(['Checking  for new files in folder ' filedir '.']);
    start(timer('TimerFcn',@timerFcn,'ExecutionMode','fixedSpacing','Period',1.0));
end

% This returns the most recently modified *bf file, recursively
% searching all subdirectories
function [fname,t] = getNewestFile(fdir)
    flist = dir(fdir);
    fname = '';
    t = 0;
    % loop through file/directory list, skip . and ..
    for k = 3 : length(flist)
        if (flist(k).isdir)
            [fn t0] = getNewestFile(fullfile(fdir,flist(k).name));
            if (t0 > t)
                fname = fn;
                t = t0;
            end
        else
            if (flist(k).datenum > t && strcmp(flist(k).name(end-1:end),'bf'))
                fname = fullfile(fdir,flist(k).name);
                t = flist(k).datenum;
            end
        end
    end
end