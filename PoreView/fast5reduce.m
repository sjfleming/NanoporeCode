function fast5reduce(filename)
%FAST5REDUCE Build reduced dataset for a fast5 Minion raw datafile

    % how many to load per segment
    nstep = 2^24;

    nchan = 512;
    numpts = 0;
    
    chans = [];
    
    % first pass, find which channels were saved
    for chan=1:nchan
        chanstr = ['/Raw/Channel_', num2str(chan), '/Signal'];
        try
            s = h5info(filename, chanstr);
        catch
            % this channel isn't present, so skip it
            continue
        end
        % save number of points, should be the same across the board
        numpts = s.Dataspace.Size;
        chans(end+1) = chan;
    end
    
    minchan = min(chans);
    maxchan = max(chans);
    nchan = maxchan-minchan+1;
    
    % now many reduced points do we aim for?
    nred = 2^19;
    % how many times to halve the data, to get at most nred points
    nhalve = ceil(log2(numpts)) - round(log2(nred));
    % how many points this leaves us with
    nred = floor(numpts/(2^nhalve));
        
    fprintf('\n\nBuilding reduced data with %d points in channels %d-%d:  0%%\n',nred,minchan,maxchan);
    
    ntotal = nchan*nred;
    
    % output mat-file
    mf = matfile([filename '_red.mat']);
    
    for chan=minchan:maxchan

        curind = 0;
        redind = 1;
        
        % output array, just this one channel
        datared = zeros(nred,1);

        while curind < numpts
            % current index and next index
            % load next data slice
            
            % are we at the end?
            thisstep = nstep;
            if (curind+nstep+1 > numpts)
                thisstep = numpts-curind-1;
            end
                
            d = fast5load(filename,[curind+1,curind+1+thisstep],chan);

            % how many points we actually got, and reduced
            np = size(d,1);
            nr = floor(np/2^nhalve);

            % check if we have too little, and pad end
            if size(d,1) < nstep
                d = [d; nan(nstep-size(d,1),size(d,2))];
            end

            % now reduce and stuff
            if nhalve > 0
                nd = 2^(nhalve+1);
                % turn into array with each column stuff and stuff
                dr = reshape(d,[nd numel(d)/nd]);
                % now each column is nd adjacent elements, alternate
                % min and max...
                %[2*numel(d)/nd 1] --> [2*numel(d)/nd/size(d,2) size(d,2)]
                d = reshape([min(dr); max(dr)],[2*numel(d)/(nd*size(d,2)), size(d,2)]);
                % seriously, just trust me guys...
            end

            % display percent loaded something something foo
            nsofar = min(ntotal, redind + (chan-minchan)*nred);
            if (np == nstep)
                fprintf('\b\b\b\b%2d%%\n',floor(100*nsofar/ntotal));
            end

            curind = curind + nstep;
            % save it into output array
            datared(redind:redind+nr-1) = d(1:nr,:);
            redind = redind+nr;
        end
        
        % now save to mat-file, just this one var
        mf.(['ch_' num2str(chan)]) = datared;
    end
    % and number of points
    mf.nred = nred;
    mf.nhalve = nhalve;
end

