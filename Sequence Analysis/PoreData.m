classdef PoreData < handle
    %POREDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Global properties
        Name        % Base name of dataset folder
        Path        % Above, including path
        NumFiles    % How many files we are looking at
        
        % Per-file parameters
        Filename        % What it sounds like (without path), for easy access
        Timestamp       % Timestamp of the file        
        Header          % ID/name of file, for later lookup
        
        Channel         % Recording channel
        Mux             % Mux setting
        ReadId          % Which read number in channel (as in filename)
        
        RawLink         % Where in fast5 file raw squiggles are stored
        RawEvents       % How many events, before splitting
        
        NumEvents       % How many events, Nx2 array (temp, comp)
        NumBases        % How many in basecalled sequence, Nx3
        NumSkips        % Number of skips from Oxford 1D summary (Nx2)
        NumStays        % Number of stays from Oxford 1D summary (Nx2)
        
        StartTime       % Of the overall event, in magic Oxford-units
        Duration        % (which are probably samples)

        ModelName       % Name of the model used in basecall
        Score           % Oxford-calculated score, Nx3
    end
    
    properties (Access=private)
        EventStore      % Internal storage of cached events/models
        EventInds       % Indices of loaded events in store (incl. negatives)
    end
    
    methods
        
        function obj = PoreData(pathname)
            
            % trim path separator from end of path
            if pathname(end) == '/' || pathname(end) == '\'
                pathname = pathname(1:end-1);
            end
                        
            if ~isempty(dir([pathname '.mat']))
                load([pathname '.mat']);
            else
                obj.Path = pathname;
                obj.AnalyzeFiles();
            end

            obj.Path = pathname;
            [~, obj.Name] = fileparts(obj.Path);
            
            obj.EventStore = [];
            obj.EventInds = [];
        end
        
        function event = getEvent(obj, ind, typ)
            
            ind = obj.getIndex(ind);
            
            filename = [obj.Path '/' obj.Filename{ind}];
            
            loc = 'template';
            if nargin == 3 && typ(1) == 'c'
                loc = 'complement';
            end
            
            try
                event = h5read(filename,['/Analyses/Basecall_2D_000/BaseCalled_' loc '/Events']);
            catch
                event = [];
                return
            end
            
            attloc = ['/Analyses/Basecall_2D_000/Summary/basecall_1d_' loc];

            shift = h5readatt(filename, attloc, 'shift');
            scale = h5readatt(filename, attloc, 'scale');
            scalesd = h5readatt(filename, attloc, 'scale_sd');
            drift = h5readatt(filename, attloc, 'drift');
            var = h5readatt(filename, attloc, 'var');
            varsd = h5readatt(filename, attloc, 'var_sd');

            model = h5read(filename,['/Analyses/Basecall_2D_000/BaseCalled_' loc '/Model']);

            model.level_mean = model.level_mean*scale + shift;
            model.level_stdv = model.level_stdv*var;
            model.sd_mean = model.sd_mean*scalesd;
            model.sd_stdv = model.sd_stdv/sqrt(varsd);
            % now drift offset to event data
            event.mean = event.mean - drift*(event.start - event.start(1));
            % and add default ref_align and ref_like
            event.ref_align = 0*event.mean;
            event.ref_like = 0*event.mean;
            % and save which one we requested
            event.type = typ;
            % and save the model too
            event.model = model;
            % and its name
            event.model.name = h5readatt(filename, attloc, 'model_file');
            while iscell(event.model.name)
                event.model.name = event.model.name{1};
            end
            % and in the model, save skip and stays
            if ~isempty(obj.NumSkips)
                if typ(1) == 't'
                    event.model.skip_prob = obj.NumSkips(ind,1)/obj.NumEvents(ind,1);
                    event.model.stay_prob = obj.NumStays(ind,1)/obj.NumEvents(ind,1);
                else
                    event.model.skip_prob = obj.NumSkips(ind,2)/obj.NumEvents(ind,2);
                    event.model.stay_prob = obj.NumStays(ind,2)/obj.NumEvents(ind,2);
                end
            else
                event.model.skip_prob = 0.04;
                event.model.stay_prob = 0.10;
            end
            event.model.extend_prob = event.model.stay_prob;
            % temp seq
            event.sequence = [];
            % flip them both for complement
            event.flipped = false;
            if typ(1) == 'c'
                event = flip_event(event);
                event.flipped = true;
            end
            % and save the (2d) sequence, while we're at it, since event
            % now points forward and whatnot
            event.sequence = obj.getSequence(ind,'2d');
            if isempty(event.sequence) && typ(1) == 't'
                event.sequence = obj.getSequence(ind,'t');
            elseif isempty(event.sequence) && typ(1) == 'c'
                event.sequence = seqrcomplement(obj.getSequence(ind,'c'));
            end
            
        end
        
        function event = getSquiggles(obj, ind)
            
            ind = obj.getIndex(ind);
            
            filename = [obj.Path '/' obj.Filename{ind}];
            
            try
                event = h5read(filename,obj.RawLink{ind});
            catch
                event = [];
                return
            end
            
        end
        
        function seq = getSequence(obj, ind, typ)
            
            seq = obj.getFastq(ind,typ);
            
            % trim out fastq quality data and leading header
            pi = find(seq=='+',1,'first');
            po = find(seq==10,1,'first');
            seq = seq(po+1:pi-2);
            
        end
        
        function seq = getFasta(obj, ind, typ)
            
            seq = obj.getFastq(ind, typ);
            % trim out fastq quality data, leave header
            pi = find(seq=='+',1,'first');
            if ~isempty(seq)
                seq = ['>' seq(2:pi-2) char(10)];
            end
            
        end
        
        function seq = getFastq(obj, ind, typ)
            
            ind = obj.getIndex(ind);
            
            filename = [obj.Path '/' obj.Filename{ind}];
            
            loc = 'template';
            if nargin == 3 && typ(1) == 'c'
                loc = 'complement';
            elseif nargin == 3 && typ(1) == '2'
                loc = '2D';
            end

            % get overall sequence data
            try
                seq = h5read(filename,['/Analyses/Basecall_2D_000/BaseCalled_' loc '/Fastq']);
            catch
                seq = [];
                return
            end
            
        end
        
        function [skips, stays] = getAlStats(obj, ind)
            
            skips = [0 0];
            stays = [0 0];
            
            ind = obj.getIndex(ind);
            % get the 2D alignment pew pew
            filename = [obj.Path '/' obj.Filename{ind}];
            ald = h5read(filename,'/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment');
            km = ald.kmer';
            % figure out where kmer does not advance
            kmstay_t = [0; all(km(1:end-1,:) == km(2:end,:),2)];
            kmstay_c = [all(km(1:end-1,:) == km(2:end,:),2); 0];
            
            % and now total stays
            stays(1) = sum(kmstay_t & (ald.template>0));
            stays(2) = sum(kmstay_c & (ald.complement>0));
            
            % and skips, a -1 where kmer does not stay
            skips(1) = sum((~kmstay_t) & (ald.template<0));
            skips(2) = sum((~kmstay_c) & (ald.complement<0));
        end
        
        function al = getAlignment(obj, ind)
            
            ind = obj.getIndex(ind);
            
            filename = [obj.Path '/' obj.Filename{ind}];
            
            try
                ald = h5read(filename,'/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment');
            catch
                al = [];
                return
            end
            al = double([ald.template, obj.NumEvents(ind,2)-ald.complement-1]);
            al(ald.complement==-1,2) = -1;
            al(al==-1) = nan;
            al = al + 1;

            al1 = al(1,1);
            for i=2:size(al,1)
                if isnan(al(i,1))
                    al(i,1) = al1;
                else
                    al1 = al(i,1);
                end
            end
            al1 = al(end,2);
            for i=size(al,1)-1:-1:1
                if isnan(al(i,2))
                    al(i,2) = al1;
                else
                    al1 = al(i,2);
                end
            end
            
        end
        
        function pores = getPores(obj)
            pores = accumarray(obj.Channel',1,[512, 1]);
        end
        
        function ind = getIndex(obj, ind)
            
            if iscell(ind)
                % cell array of names
                inds = zeros(numel(ind),1);
                for i=1:numel(ind)
                    inds(i) = find(strcmp(obj.Header, ind{i}));
                end
                ind = inds;
            elseif ischar(ind)
                % a single name
                ind = find(strcmp(obj.Header, ind));
            end
            
        end
        
        function events = getEvents(obj, evinds)
            % Update the event store with events and return
            
            keepev = zeros(size(obj.EventStore));
            
            % figure out which ones we can keep
            for i=1:numel(obj.EventStore)
                if any(obj.EventInds(i) == evinds)
                    keepev(i) = 1;
                end
            end
            % and select only those
            obj.EventInds = obj.EventInds(keepev > 0);
            obj.EventStore = obj.EventStore(keepev > 0);
            % now go the other way, and add events
            for i=1:numel(evinds)
                % check if we covered it already
                if any(obj.EventInds == evinds(i))
                    continue
                end
                % and if not, load and add
                event = obj.getEvent(abs(evinds(i)),'t');
                if ~isempty(event)
                    % flip if it's reverse-aligned
                    % this flips internal sequence too
                    if evinds(i) < 0
                        event = flip_event(event);
                    end

                    % now save it
                    if isempty(obj.EventStore)
                        obj.EventStore = event;
                        obj.EventInds = evinds(i);
                    else
                        obj.EventStore(end+1) = event;
                        obj.EventInds(end+1) = evinds(i);
                    end
                end
                
                % do the same for the complement
                event = obj.getEvent(abs(evinds(i)),'c');
                if ~isempty(event)
                    % flip if it's reverse-aligned
                    if evinds(i) < 0
                        event = flip_event(event);
                    end

                    % now save it
                    if isempty(obj.EventStore)
                        obj.EventStore = event;
                        obj.EventInds = evinds(i);
                    else
                        obj.EventStore(end+1) = event;
                        obj.EventInds(end+1) = evinds(i);
                    end
                end
            end
            %{
            if nargin > 2
                fprintf('Using existing alignments\n');
                for i=1:numel(obj.EventInds)
                    evind = obj.EventInds(i);
                    al = als(evinds == evind,:);
                    evn = numel(obj.EventStore(i).ref_align);
                    al_m = (al(2)-al(1))/(al(4)-al(3))/evn;
                    al_b = al(1) - al_m * al(3);
                    obj.EventStore(i).ref_align = round(al_m * (1:evn)' + al_b);
                    obj.EventStore(i).ref_align = min(obj.EventStore(i).ref_align, al(2));
                    obj.EventStore(i).ref_align = max(obj.EventStore(i).ref_align, al(1));
                end
            end
            %}
            % and return
            events = obj.EventStore;
            
        end
        
        function seqs = getSequences(obj, evinds, typ)
            % Returns a cell array of sequences for specified event indices
            % of class typ (if available)
            
            seqs = {};
            
            for i=1:numel(evinds)
                s = obj.getSequence(abs(evinds(i)),typ);
                if isempty(s)
                    continue
                end
                % reverse only if not complement/reversed thing
                if xor(typ(1) == 'c', evinds(i) < 0)
                    s = seqrcomplement(s);
                end
                seqs{end+1} = s;
            end
        end
        
    end
    
    methods (Hidden)
        
        function AnalyzeFiles(obj)
            
            fns = dir([obj.Path '/*.fast5']);

            if isempty(fns)
                error('No files found!')
            end
            
            obj.NumFiles = numel(fns);

            tic;
            
            for i=1:numel(fns)
                
                filename = [obj.Path '/' fns(i).name];
                obj.Filename{i} = fns(i).name;
                obj.Timestamp(i) = fns(i).datenum;

                try
                    loc = '/Analyses/Basecall_2D_000/Configuration/general';
                    obj.Header{i} = h5readatt(filename,loc,'tag');
                    while iscell(obj.Header{i})
                        obj.Header{i} = obj.Header{i}{1};
                    end
                    obj.Channel(i) = str2double(h5readatt(filename,loc,'channel'));
                    obj.ReadId(i) = str2double(h5readatt(filename,loc,'read_id'));
                    
                    
                    hi = h5info(filename,'/Analyses/Basecall_2D_000');
                    obj.RawLink{i} = ['/' hi.Links.Value{1}];
                    hd = h5info(filename,obj.RawLink{i});
                    obj.RawEvents(i) = hd.Dataspace.Size;
                    attlink = obj.RawLink{i}(1:find(obj.RawLink{i}=='/',1,'last')-1);
                    obj.StartTime(i) = h5readatt(filename,attlink,'start_time');
                    obj.Duration(i) = h5readatt(filename,attlink,'duration');
                    
                    obj.Mux(i) = str2double(h5readatt(filename,attlink,'start_mux'));
                catch
                    fprintf('Error in file %s\n',filename);
                    continue
                end

                locs = {'template','complement'};
                for j=1:2
                    try
                        loc = ['/Analyses/Basecall_2D_000/Summary/basecall_1d_' locs{j}];
                        obj.ModelName{i,j} = h5readatt(filename,loc,'model_file');
                        obj.NumBases(i,j) = h5readatt(filename,loc,'sequence_length');
                        obj.Score(i,j) = h5readatt(filename,loc,'mean_qscore');
                        obj.NumSkips(i,j) = h5readatt(filename,loc,'num_skips');
                        obj.NumStays(i,j) = h5readatt(filename,loc,'num_stays');
                        
                        hd = h5info(filename,['/Analyses/Basecall_2D_000/BaseCalled_' locs{j} '/Events']);
                        obj.NumEvents(i,j) = hd.Dataspace.Size;
                    end
                end
                try
                    obj.NumBases(i,3) = h5readatt(filename,'/Analyses/Basecall_2D_000/Summary/basecall_2d','sequence_length');
                    obj.Score(i,3) = h5readatt(filename,'/Analyses/Basecall_2D_000/Summary/basecall_2d','mean_qscore');
                end

                if mod(i,100) == 0
                    t = toc;
                    trem = (numel(fns)-i)*t/i;
                    fprintf(1,'%0.2f minutes remaining...\n',trem/60);
                end
            end
            
            save([obj.Path '.mat'],'obj');
        end
        
    end
    
end

