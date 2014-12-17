classdef SignalData < handle
    %SIGNALDATA: Class wrapper for streaming signals (specifically abfs)
    %
    % ---------------------------------------------------------------------
    %
    % SignalData Methods:
    %   SignalData(fname) - Initialize class on a file, if it exists
    %   getViewData(trange) - Return reduced or full data in a time range
    %   get(inds,sigs) - Return full data in specified index range
    %   getByTime(t0,t1) - Return full data in specified time range
    %   addVirtualSignal(fun,name,srcs) - Add a virtual signal function
    %   getSignalList() - Get names of all accessible signals
    %   findNext(fun,istart) - Find next instance of logical 1
    %   findPrev(fun,istart) - Find previous instance of logical 1
    %
    % ---------------------------------------------------------------------
    %
    % SignalData Properties:
    %   filename - Loaded name of file
    %   ndata - Number of data points (per signal)
    %   nsigs - Number of signals in file
    %   si - Sampling interval, seconds
    %   tstart - Start time of file, set to 0
    %   tend - End time of file
    %   header - Header struct from abf file
    %
    % ---------------------------------------------------------------------
    %
    % About Data:
    %   All data, always, anywhere, returned or used by SignalData, always
    %   has columns [time, sig1, sig2, ...], where the signals represent
    %   the different pieces of data present in the file (or added via
    %   virtual signals, see below). So, column/signal 1 is always Time.
    %   This is useful for backing out the absolute time/index of data that
    %   you are processing in a small chunk:
    %   
    %       d = obj.get(5000:5050,:);       % same as obj.get(5000:5050)
    %       ind = do_some_processing(d);    % ind is 1...50 since d has 50 points
    %       t = d(ind,1);                   % t is the absolute time of event in ind
    %                                
    % About Reduced Data:
    %   The reduced data calculated by SignalData subsamples the entire
    %   file to generate a reduced version consisting of 500k points. The
    %   subsampled (reduced) data is a series of points that alternate
    %   between the min and max value of the points they replace. In other
    %   words, reduced pt. 1 is eg. min(1:1000) and pt. 2 is max(1:1000).
    %   This way, features such as events aren't lost through subsampling.
    %   The reason min and max are alternated is so that the full range is
    %   visible when subsampled data is plotted (which appears in a lighter
    %   color). The reduced data is saved to a file, if possible, so that
    %   the next time you load a data file, the reduced data appears right
    %   away.
    % 
    % About Caching:
    %   Accessing the data using obj.get and obj.getByTime give you the
    %   full (non-reduced) data in a given range. SignalData keeps a
    %   certain amount of data (1 million points or so) in memory (called a
    %   'cache'), so if you access the points 5000:5050 and then 5050:5100,
    %   it will not read any data from disk the second time, only when you
    %   request points outside of what is loaded into memory (say, 1e6:1e6+50). 
    %   This way, you can access the file as if it were all loaded, without
    %   having to worry about accessing the disk every time.
    %   You can request as many points as you like, if you want to load
    %   more than the cache normally holds - it is up to you to ensure that 
    %   you don't request the entire file by accident.
    %   
    % About Virtual Signals: 
    %   A virtual signal is a filter you have written that acts on data
    %   with columns [time, sig1, sig2, ...], returing an array of the same
    %   form but processed by the filter. The virtual signal is applied to
    %   the full data whenever it is loaded in from the cache, and appears
    %   as a new signal column whenever you access the data (eg. obj.get()).
    %   Basically, a virtual signal is a filter that will appear as if you
    %   have a filtered version of the entire file. What function you apply
    %   is up to you, whether it is a high/low/bandpass filter, or median
    %   filter, or one that replaces a time range of points with their
    %   average.
    %   One important note is that the program attempts to apply the filters
    %   to the reduced data as well, which can give funny results depending
    %   on how well the filters are suited to taking such data (for
    %   example, highpass and lowpass work fine, but median filters make
    %   short events completely disappear).
    %
    % ---------------------------------------------------------------------
    %
    %
    %
    
    % make it so these don't get screwed up
    properties (SetAccess=immutable)
        filename    % filename we are working with
        ndata       % number of points
        nsigs       % number of signals (not including time)
        si          % sampling interval
        tstart      % start time of file, set to 0
        tend        % end time of file
        header      % original abf info header
    end
    
    % can't change these from the outside
    properties (SetAccess=private, Hidden=true)
        nred        % number of points to store in the reduced array
        datared     % subsampled data, in min-max form
                    % gets updated as virtual stuff changes
        
        cstart      % starting point of loaded cached data
        cend        % endpoint of loaded cached data
        dcache      % cached data that we're working with
        
        nvsigs      % how many virtual signals? this is the total number (incl. multiple outputs)
                    % such that data has 1 + nsigs + nvsigs columns
                    
        vnames      % cell array names of virtual signals
        vfuns       % cell array functions for virtual signals
        vsrcs       % cell array, which columns get processed by each one
        vdsts       % cell array, which columns get written to
    end
    
    methods
        function obj = SignalData(fname)
            % obj = SignalData(filename) - Creates class based on specified file
            %   Builds reduced data set if it doesn't yet exist, and tries
            %   to save it. If the file isn't loaded properly, resulting
            %   obj.ndata is set to -1, or -2 if it's an IV curve.
            
            % start working!
            obj.filename = fname;
            
            % try to load file, see if we got it right
            try
                % first, load some info on the file               
                if obj.filename(end-2) == 'a'
                    [~,~,h]=abfload(obj.filename,'info');
                else
                    disp(['Loading cbf file ' fname '...'])
                    [~,h]=cbfload(obj.filename,'info');
                end
            catch
                fprintf(2,'Failed to load file %s!\n',obj.filename);
                obj.ndata = -1;
                return
            end
            
            % check if it's an IV curve, and cry if it is
            try
                % abf version
                if isfield(h,'lSynchArraySize') && h.lSynchArraySize > 0
                    obj.ndata = -2;
                    return
                end
                % cbf version
                if isfield(h,'type') && ~strcmp(h.type,'Continuous')
                    obj.ndata = -2;
                    return
                end
            catch
                fprintf(2,'Unrecognized file!\n');
                return
            end
            
            % clear virtual signal parts
            obj.nvsigs = 0;
            obj.vnames = {};
            obj.vfuns = {};
            obj.vsrcs = {};
            
            obj.header = h;
            obj.nred = 0;
            
            if ~isfield(h,'type')
                % abf version
                
                obj.si = h.si*1e-6;
                % knock a couple points off the end, just to prevent
                % bizarre off-by-one errors...?
                obj.ndata = h.dataPtsPerChan - 2;
                obj.tstart = 0; % dunno how to get actual start from abf
                obj.tend = obj.si*(obj.ndata-1);
                obj.nsigs = h.nADCNumChannels;
                % which signals to reduce
                redsrc = 1:(1+obj.nsigs);
            else
                % cbf version
                
                obj.si = h.si;
                obj.ndata = h.numPts - 2;
                obj.tstart = 0;
                obj.tend = obj.si*(obj.ndata-1);
                obj.nsigs = h.numChan;
                
                % add a default virtual signal for host low-pass filter
                f = h.hostFilterHz;
                obj.addVirtualSignal(@(d) filt_lp(d,4,f),...
                    sprintf('Low-Pass %d',f),[2 3]);
                % use low-pass filtered for reduced data...
                redsrc = [1 6 7 4 5];
            end
            
            % set cache to default values
            obj.cstart = 0;
            obj.cend = 0;
            obj.dcache = [];

            % try to load subsampled file
            try
                tmp = load([obj.filename  '_red.mat'],'red');
                obj.datared = tmp.red;
                obj.nred = size(obj.datared,1);
                fprintf('\nLoaded reduced data from %s_red.mat.\n',obj.filename);
            catch
                % CHANGE THIS LINE TO CHANGE HOW MANY REDUCED POINTS WE
                % HAVE WHEEEEEE
                obj.nred = 5e5;
                
                % check if we have few enough points to not need reduced
                % data, which is just the above value times a constant
                if (obj.ndata < 4*obj.nred)
                    % the rest of the program will know that this means
                    % that there is no reduced data being used
                    obj.nred = 0;
                    obj.datared = [];
                    fprintf('\nNo reduced dataset needed.\n');
                    return;
                end
                
                fprintf('\n\nBuilding reduced data with %d points -  0%%\n',obj.nred);
                obj.datared = zeros(obj.nred,obj.nsigs+1);
 
                % go through entire file and build it
                % this function returns index into full data of ith point
                fullIndex = @(ind) floor((obj.ndata-1)*ind/obj.nred);
                obj.datared(:,1) = obj.si*fullIndex(1:obj.nred);
                % full points per reduced point
                dind = floor(obj.ndata/obj.nred);
                % number of red. points to do per file loading step
                % load ~ 4 sec at a time, should be an even number
                nstepred = 2*floor(1/obj.si/dind);
                
                ired = 0;
                
                while (ired < obj.nred)
                    % current index and next index
                    inext = ired + nstepred;
                    if (inext > obj.nred)
                        inext = obj.nred;
                    end
                    
                    % load next data slice
                    d = obj.getData(fullIndex(ired),fullIndex(inext));
                    % and use only the parts we said we would
                    d = d(:,redsrc);
                    % generate array of indices - each el corresponds to
                    % reduced point index (yeah this is a bit hacky oh well)
                    % eg. this array ends up as [1 1 1 1 1 2 2 2 2 2 3 etc]
                    inds = floor(linspace(1,0.999+(inext-ired)*0.5,size(d,1)))';
                    
                    % get the min and max values for each signal
                    for i=2:obj.nsigs+1
                        mins = accumarray(inds,d(:,i),[nstepred/2,1],@min);
                        maxs = accumarray(inds,d(:,i),[nstepred/2,1],@max);
                        np = 2*size(mins,1);
                        % make alternating min-max data
                        % and write it to output array
                        obj.datared(ired+1:ired+np,i) = reshape([mins maxs]',[np 1]);
                    end
                                        
                    ired = ired + nstepred;
                    
                    % display percent loaded something something foo
                    idisp = min(obj.nred,ired);
                    fprintf('\b\b\b\b%2d%%\n',floor(100*idisp/obj.nred));
                end
                
                % trim off extra points
                obj.datared = obj.datared(1:obj.nred,:);
                
                % and save the data
                try
                    red = obj.datared;
                    save([obj.filename  '_red.mat'],'red');
                    fprintf('\nDone, saved to %s_red.mat.\n',obj.filename);
                catch
                    fprintf(2,'\nCould not save reduced data to %s_red.mat!\n',obj.filename);
                end
                
                % make sure everything is tidy
                obj.updateVirtualData(true);
            end
        end

        
        function [d, isred] = getViewData(obj,trange)
            % [data, isReduced] = obj.getViewData([tstart tend])
            %   Returns reduced or full data in a specified time range,
            %   with the full one being returned once it wouldn't kill the
            %   computer. Also tells you if it's using reduced or full
            %   version.
            
            dt = trange(2)-trange(1);
            
            % trim the time range if it's too big
            if (trange(1) < 0)
               trange(1) = 0;
            end
            if (trange(2) < 0)
               trange(2) = 0;
            end
            if (trange(1) > obj.tend)
                trange(1) = obj.data.tend;
            end
            if (trange(2) > obj.tend)
                trange(2) = obj.tend;
            end
            
            % reduced sampling interval
            redsi = (obj.tend-obj.tstart)/obj.nred;
            % number of points from reduced set we'd use
            nr = dt/redsi;
            % number of points from full set we would be using
            nfull = dt/obj.si;
            % only use reduced one if it wouldn't be visible (nr>1500)
            % and if there would be too many regular points (nfull>nred)
            if nfull > obj.nred && nr > 1500
                % use reduced
                inds = floor(trange/redsi);
                % already contains virtual data
                d = obj.datared(inds(1)+1:inds(2),:);
                % and set our using reduced flag
                isred = true;
            else
                % use full, if we have virtual signals the array will
                % be bigger and stuff and junk
                pts = floor(trange/obj.si);
                d = obj.getData(pts(1),pts(2));
                % and set flag
                isred = false;
            end
        end
    
        
        function d = get(obj,pts,sigs)
            % data = obj.get(inds) - Get all signals in range inds
            % data = obj.get(inds,sigs) - Get specified signals in range inds
            %   Returns data by index. obj.get(5:43), obj.get(5:43, 1:4)
            %   etc. Returns in the full specified range, so for example
            %   obj.get(5:43) is the same as obj.get( [5,43] ).
            %   If points outside data file limit are requested, will
            %   trim and possibly return zero points.
        
            % if we didn't specify which signals, take all of them
            if nargin < 3
                sigs = ':';
            end
            % get the data, including time
            d = obj.getData(min(pts),max(pts));
            % return only requested signals
            d = d(:,sigs);
        end
        
        
        function d = getByTime(obj,t0,t1)
            % data = obj.getByTime(t0, t1)
            % data = obj.getByTime([t0 t1])
            %   Return data points in specified time range, if possible.
            
            if nargin == 3
                t0 = [t0 t1];
            end
            pts = floor(t0/obj.si);
            d = obj.getData(min(pts),max(pts));
        end
        
        
        function dst = addVirtualSignal(obj, fun, name, src)
            % dst = obj.addVirtualSignal(fun) - Add a function as a virtual signal
            % dst = obj.addVirtualSignal(fun, name) - Give it a name, too
            % dst = obj.addVirtualSignal(fun, name, src) - And which signals to pass to the function
            %   If signal with name exists, it gets replaced. Either way, 
            %   returns which columns the virtual signal appears as.
        
            % if we didn't specify source, do time + orig. signals
            if (nargin < 4)
                src = 1:obj.nsigs+1;
            else
                % and if we did, make sure 1 is on there
                if isempty(find(src==1,1))
                    src = [1 src];
                end
            end
            
            % how many new signals are we adding?
            nadd = length(src) - 1;
            
            % check if this one exists already, if we were given a name
            if nargin > 2 && any(ismember(obj.vnames, name))
                i = find(ismember(obj.vnames, name),1);
                % if it exists, we keep the destination the same
                dst = obj.vdsts{i};
            else
                % or add a new one
                i = length(obj.vnames) + 1;
                obj.nvsigs = obj.nvsigs + nadd;
                % dst is the last nadd columns
                dst = (-nadd+1:0) + obj.nvsigs + obj.nsigs + 1;
            end
            
            % if no name given, make one up
            if (nargin < 3)
                name = sprintf('Virtual %d',i);
            end
            
            obj.vnames{i} = name;
            obj.vfuns{i} = fun;
            obj.vsrcs{i} = src;
            obj.vdsts{i} = dst;
            
            % and now that we've added it, make sure reduced and cached data's good
            obj.updateVirtualData(true);
        end
               
        
        function siglist = getSignalList(obj)
            % names = obj.getSignalList()
            %   Return a list of accessible signals, in order they appear,
            %   as a cell array.
            
            % first signal is always time
            siglist = {'Time'};
            for i=1:obj.nsigs
                if isfield(obj.header,'recChNames')
                    siglist{end+1} = obj.header.recChNames{i};
                else
                    siglist{end+1} = obj.header.chNames{i};
                end
            end
            % for virtual signals, append the filter name
            for i=1:length(obj.vnames)
                src = obj.vsrcs{i};
                for j=src(2:end)
                    siglist{end+1} = sprintf('%s (%s)',obj.vnames{i},siglist{j});
                end
            end
        end
        
        
        function ind = findNext(obj,fun,istart)
            % ind = obj.findNext(fun)
            % ind = obj.findNext(fun, istart)
            %   Finds next instance of logical 1, starting at index, if
            %   specified.
            
            % we don't need to specify istart
            if (nargin < 3)
                istart = 0;
            end
            
            % number of points to step by, hard-coded for now
            maxPts = 1e5;
            
            % loop and find next index of a logical 1
            while 1
                d = obj.get(istart:istart+maxPts);
                
                % check if we have hit the end of the file?
                if isempty(d)
                    % then give up and cry
                    ind = -1;
                    return
                end
                
                % find the index! (if we have one)
                ind = find(fun(d),1,'first');
                
                % did we find a logical 1?
                if ~isempty(ind)
                    % shift index and return it
                    ind = ind + istart - 1;
                    return
                end
                
                istart = istart + maxPts;
            end
        end
        
        function ind = findPrev(obj,fun,istart)
            % ind = obj.findPrev(fun)
            % ind = obj.findPrev(fun, istart)
            %   Finds previous instance of logical 1, starting at index, if
            %   specified.
            
            % we don't need to specify istart
            if (nargin < 3)
                istart = obj.ndata;
            end
            
            % number of points to step by, hard-coded for now
            maxPts = 1e5;
            
            % loop and find next index of a logical 1
            while 1
                i0 = istart-maxPts;
                d = obj.get(i0:istart);
                
                % check if we have hit the end of the file?
                if isempty(d)
                    % then give up and cry
                    ind = -1;
                    return
                end
                
                % find the index! (if we have one)
                ind = find(fun(d),1,'last');
                
                % did we find a logical 1?
                if ~isempty(ind)
                    % shift index and return it
                    ind = ind + i0 - 1;
                    return
                end
                
                istart = istart - maxPts;
            end
        end
        
    end
        
    
    % internal functions go here
    methods (Access=private, Hidden=true)
        
        function d = getData(obj, ptstart, ptend)
            % data = obj.getData(ptstart,ptend)
            %   Returns data in the specified point range, from cache. Also
            %   updates the cache if necessary. This is for internal use.
        
            % check bounds first thing
            if (ptstart < 0 || ptend > obj.ndata-1 || ptend<ptstart)
                fprintf(2,'Invalid points %d:%d requested\n',int64(ptstart),int64(ptend));
            end
            ptstart = max(0,ptstart);
            % keep one point away from end, cause this is actually one too
            % many points, and we can't really request it
            ptend = min(obj.ndata-1,ptend);
            
            if (ptstart < obj.cstart || ptend >= obj.cend)
                % cache miss, load a new cache
                % size range requested
                dpt = (ptend-ptstart);
                % conservatively load a million points (or more if needed)
                if (dpt < 1e6)
                    % extend loading range to ~1 million
                    obj.cstart = round(ptstart - (1e6-dpt)/2);
                    obj.cend = round(ptend + (1e6-dpt)/2);
                else
                    % or just by 10 each way to avoid indexing errors etc
                    obj.cstart = ptstart-10;
                    obj.cend = ptend+10;
                end
                % and check bounds
                obj.cstart = max(obj.cstart,0);
                obj.cend = min(obj.cend,obj.ndata);
                
                % load file, add in a 'cheat point' at the end just to make
                % sure we get everything
                if ~isfield(obj.header,'type')
                    % abf version!
                    d = abfload(obj.filename,'start',obj.cstart*obj.si,'stop',...
                        (obj.cend+1)*obj.si,'verbose',0);
                else
                    % cbf version
                    d = cbfload(obj.filename,[obj.cstart,(obj.cend+1)]);
                end
                %fprintf('Loaded %d points (%d-%d) into the cache\n   ',size(obj.dcache,1),floor(obj.cstart),floor(obj.cend));
                
                % make empty cache
                obj.dcache = zeros(size(d,1),1 + obj.nsigs + obj.nvsigs);
                
                % add time data on to cache, as first col
                npts = size(obj.dcache,1);
                ts = obj.si*((1:npts)-1+obj.cstart);
                % set, along with loaded data
                obj.dcache(:,1:obj.nsigs+1) = [ts' d];
                
                % and update the virtual signals, but not for reduced
                obj.updateVirtualData(false);
            end
            % now we definitely have the points
            % +1 for Matlab's 1-indexed arrays ugh
            pts = int64(ptstart - obj.cstart+1);
            pte = int64(ptend - obj.cstart+1);
            
            d = obj.dcache(pts:pte,:);
        end
        
        function updateVirtualData(obj, dored)
            % obj.updateVirtualData(dored)
            %   Is exactly what it sounds like. Updates all of the internal
            %   virtual data, in full and reduced, if requested.
            
            if (dored && obj.nred > 0)
                % set original reduced data aside
                d = obj.datared(:,1:obj.nsigs+1);
                % make the new one
                obj.datared = zeros(obj.nred, 1+obj.nsigs+obj.nvsigs);
                % set the originals
                obj.datared(:,1:obj.nsigs+1) = d;
            end

            % do the same with the cache
            if ~isempty(obj.dcache)
                d = obj.dcache(:,1:obj.nsigs+1);
                obj.dcache = zeros(size(obj.dcache,1), 1+obj.nsigs+obj.nvsigs);
                obj.dcache(:,1:obj.nsigs+1) = d;
            end

            % and apply virtual signal functions to cached and reduced data
            for i=1:length(obj.vnames)
                % which columns to write to?
                dst = obj.vdsts{i};
                % and which to read from
                src = obj.vsrcs{i};
                fun = obj.vfuns{i};
                % and execute it
                if ~isempty(obj.dcache)
                    % just a check to make sure we get the right number
                    % of columns from the virtual functions
                    A = fun(obj.dcache(:,src));
                    obj.dcache(:,dst) = A(:,(end-length(dst)+1):end);
                end
                if (dored && obj.nred > 0)
                    A = fun(obj.datared(:,src));
                    obj.datared(:,dst) = A(:,(end-length(dst)+1):end);
                end
            end
        end
    end
end

