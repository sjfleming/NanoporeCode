function DNAevent = find_events(pv)
%FIND_EVENTS Finds events and returns an array of DNAevent structs
%   events = find_events(pv)
%   This is just an example of how to use PoreView programmatically.

    % Alternatively, instead of passing a PoreView instance to the code,
    % you could just initialize one here:
    % pv = PoreView(filename);
    % and then add any virtual signals you need, set up the panels, and
    % then start finding events on your merry way...

    DNAevent = [];

    thresh = 0.07;
    
    % use the signal specified in the second panel of PoreView for the
    % event finding. you'd probably just hardcode this in your own code.
    sig = pv.psigs(2).sigs;

    % loop through entire file, a bit at a time
    curind = 0;
    
    while 1
        % find next data exceeding threshold, stepping current index
        curind = pv.data.findNext(@(d) d(:,sig) > thresh, curind);
        
        % if we didn't find any, we're done with the file
        if curind < 0
            break
        end
        
        imin = curind;
        % find the end of the event
        imax = pv.data.findNext(@(d) d(:,sig) < 0.75*thresh,curind);
        
        % make sure we have an end for the event
        if imax < 0
            break
        end
        
        % shift event by one sample in each directon to get whole event
        imin = imin-1;
        imax = imax+1;
        
        % and move curind too
        curind = imax;
        
        % event? maybe event?
        ts = pv.data.si*[imin imax];
        
        % time range to view, extended past the event a bit
        viewt = [ts(1)-0.001, ts(2)+0.001];
        
        % we've decided we have a possibly good event, then
        dna = [];
        
        % store the data we want, including times
        % note that we're only grabbing the signal we're analyzing
        dna.data = pv.data.get(imin:imax,[1 sig]);

        % and the start and end times for the event
        dna.tstart = ts(1);
        dna.tend = ts(2);

        % the average current blockage
        dna.blockage = abs(mean(dna.data(:,2)));
        
        % now query on-screen, see what we think
        % first, zoom in
        pv.setView(viewt);
        % then, draw some stuff
        h = pv.getAxes(2);
        plot(h, [viewt(1) ts(1) ts(1) ts(2) ts(2) viewt(2)],...
            [0 0 dna.blockage dna.blockage 0 0],'r');
        % this line ignores the stuff you drew, in case you're wondering
        pv.autoscaleY();
        % do this just for kicks
        pv.setCursors(ts);
        
        % if you want it to run automatically, this would be the part you
        % would change, or something
        k = pv.waitKey();
        % clear, and force a redraw. this way, you don't accidentally press
        % keys twice because you don't know if it's thinking or not.
        pv.clearAxes();
        pause(0.01);
        
        % handle key input
        if (k == 'q')
            return
        elseif (k ~= 'y')
            continue
        end
        
        % stupid matlab structs bug
        if isempty(DNAevent)
            DNAevent = dna;
        else
            DNAevent(end+1) = dna;
        end
    end
end

