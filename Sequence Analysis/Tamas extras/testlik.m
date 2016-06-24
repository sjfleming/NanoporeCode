function str = testlik(  )

    
    mex align_like.cpp
    mex viterbikd.cpp
 

    pathname = 'C:\Minion\Lambda-freezein';
    
    filename = [pathname '_template.fasta.bam'];

    % first, load bam stuff into memory
    try
        h=baminfo(filename,'ScanDictionary',true);
    catch
        fprintf(2,'Not found: %s\n',filename);
    end
    bm = BioMap(filename, 'SelectReference', h.ScannedDictionary{1});

    
    lambda = fastaread('.\References\Lambda_NEB.fasta');
    lambda = lambda.Sequence;
    
    pd = PoreData(pathname);

    % get headers that align with the file
    refstart = 10000;
    refend = 11000;
    bminds = bm.getIndex(refstart,refend,'overlap',refend-refstart);
    headers = bm.Header(bminds);
    % and get their event indices
    evinds = zeros(size(headers));
    for i=1:numel(headers)
        evinds(i) = pd.getIndex(headers{i});
    end
    
    fprintf('Starting with %d strands...\n',numel(headers));

    % loop through and add events, t and c
    for i=1:numel(evinds)
        event = pd.getEvent(evinds(i),'t');
        event.ref_align = 0*event.mean;
        
        if bm.Flag(bminds(i)) > 0
            event = flip_event(event);
        end
        
        % now save it
        events(i) = event;
    end
    
    % patchy code to make deviations not matter
    for i=1:0%numel(events)
        foo = 0*events(i).model.sd_stdv;
        events(i).model.sd_stdv = foo + 1/sqrt(2*pi);
        events(i).model.sd_mean = foo + 0;
        events(i).stdv = 0*events(i).stdv;
    end
    

    lstr = lambda(refstart:refend);
    % pick one of the strands to start with
    strind = evinds(bm.Flag(bminds) == 0);
    str = pd.getSequence(strind(end-2),'2d');
    [~,p,st] = swalign(str,lstr,'Alphabet','NT');
    % how many from first strand in this alignment?
    ns = sum(p(1,:) ~= '-');
    %ns = 1000;
    fprintf('Using strand from %d to %d\n',st(1),st(1)+ns);
    str = str(st(1):st(1)+ns);
    
    %str = lstr;
    
    str_orig = str;
    str = nt2int(str);

    tic
    
    laststr = str;
    
    states = get_states(str);
    
    
    rng(13);
    
    allevents = events(randperm(numel(events)));
    
    events = allevents(1:end);
    
    tic
    strind = 1;
    
    score = -1e100;
    
    function better = runstrands(curstr)
        curscore = 0;
        for k=1:numel(events)
            [s, events(k)] = align_like(events(k), get_states(curstr));
            curscore = curscore + s;
        end
        better = false;
        if curscore > score
            str = curstr;
            score = curscore;
            better = true;
        end
    end

    numloops = 5;
    
    while 1
        % test mutations
        for i=1:4
            for j=1:4
                newstr = str;
                newstr(strind:strind+1) = [i j];
                if runstrands(newstr)
                    strind = max(strind-3,1);
                end
            end
        end
        % and deletions
        newstr = str;
        newstr = newstr(1:numel(newstr) ~= strind);
        if runstrands(newstr)
            strind = max(strind-3,1);
        end
        
        % and insertions
        
        for i=1:-1
            newstr = str;
            newstr = [newstr(1:strind) i newstr(strind+1:end)];
            if runstrands(newstr)
                strind = max(strind-3,1);
            end
        end
        
        fprintf('%d: %0.1f | %f | %d\n',strind,score,swalign(lstr,int2nt(str),'Alphabet','NT'), ...
            numel(str));
        strind = strind + 1;
        
        if (strind >= numel(str) - 1)
            numloops = numloops - 1;
            strind = 1;
            swalign(lambda,str,'Alphabet','NT')
            if (numloops == 0)
                break
            end
        end
    end
    toc

    str = int2nt(str);
    
    [s,p] = swalign(lambda,str,'Alphabet','NT')
    
    [s,p] = swalign(lambda,str_orig,'Alphabet','NT')

end

