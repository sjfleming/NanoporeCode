function testparams( )

    makeall    

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

    % loop through and add events
    for i=1:numel(evinds)
        [event, model] = pd.getEvent(evinds(i),'t');
        
        if bm.Flag(bminds(i)) > 0
            event = flip_event(event);
            model = flip_model(model);
        end
        
        % now save it
        events(i) = event;
        models(i) = model;
    end
    
    
    fprintf('Now doing compliments....\n');
    
    for i=1:numel(evinds)
        [event, model] = pd.getEvent(evinds(i),'c');
        
        if isempty(event)
            continue
        end
        
        % flip in opposite situation
        if bm.Flag(bminds(i)) == 0
            event = flip_event(event);
            model = flip_model(model);
        end
        
        % now save it
        events(end+1) = event;
        models(end+1) = model;
    end
    
    % patchy code to make deviations not matter
    for i=1:numel(events)
        foo = 0*models(i).sd_stdv;
        models(i).sd_stdv = foo + 1/sqrt(2*pi);
        models(i).sd_mean = foo + 0;
        events(i).stdv = 0*events(i).stdv;
        % and let's put a lower bound on the level deviations
        %lmean = mean(models(i).level_stdv);
        %models(i).level_stdv = max(models(i).level_stdv,lmean);
    end
    
    % use all of the events for this bit
    % and start with a pretty good sequence
    load('bestseq.mat')

    lstr = lambda(refstart:refend);
    [~,p] = swalign(seq,lstr,'Alphabet','NT');
    seq = p(1,:);
    seq = seq(seq ~= '-');
    
    % get rid of homopolymer repeats, for now
    bases = 'ACGT';
    for i=1:4
        strshort = repmat(bases(i),[1,5]);
        lstr = regexprep(lstr,[strshort '+'],strshort);
    end
        
    % launch a viewer for looking at the sequence as we mutate it
    pv = PoreView();
    pv.setView([0.5 1000.5]);
    zoop = onCleanup(@() delete(pv.fig));
    
    alparams = [];
    alparams.stripe_width = 250;
    alparams.comp_cutoff = 50;
    alparams.insert_prob = 0.001;
    alparams.stay_prob = 0.05;
    alparams.extend_prob = 0.08;
    alparams.skip_prob = 0.10;
    alparams.lik_offset = 4.5;
    
    % seed alignment
    [~,~,events] = align_likes(events, models, seq, alparams);
    
    for i=1:1
        % plot the stats as we go
        stats = align_stats(seq, events, models, lstr, alparams);
        plot_stats(seq, stats, events, models, pv);
        plot_events(events,models,seqtostates(seq),845);
        pause(0.1);
        
        % and measure scores
        seqscore = sum(align_likes(events, models, seq, alparams));
        lscore = sum(align_likes(events, models, lstr, alparams));
        
        fprintf('Seq/Lambda: %0.1f/%0.1f\n',seqscore,lscore);
    end
end

