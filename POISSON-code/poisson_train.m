function poisson_train()

    make all
    
    % load the parameters file, if one exists
    try
        params = PoreParams.Load('CS_params.conf');
    catch
        params = PoreParams.Default();
    end
    
    % and the DNA CS sequence
    refseq = fastaread('References/DNA_CS.fasta');
    refseq = refseq.Sequence;
    
    parc = parcluster('local');
    if isunix()
        parc.NumWorkers = 16;
    else
        parc.NumWorkers = 4;
    end
    nt = parc.NumWorkers;
    
    try
        parpool(parc);
    end
    
    % now let's load the data and get us some events/sequence
    if ~isunix()
        pd = PoreData('C:\Minion\DNA_CS1');
    else
        pd = PoreData('~/Minion/DNA_CS1');
    end
    % keep only those where template, complement, and 2d have enough bases
    evinds = find(min(pd.NumBases,[],2)>3000);
    
    % how many strands are we using for the training?
    coverage = 10;
    % take a subset from all aligned headers
    events = pd.getEvents(randsubset(evinds,coverage));
    % order the events
    events = order_events(refseq,events);
    % seed the alignments 
    events = seedaligns(refseq,events,params);
    events = seedaligns(refseq,events,params);
    events = seedaligns(refseq,events,params);
    
    % now using default/current params, take refseq down to baseline
    seq = refseq;
    for i=1:3
        [seq,events] = MutateSequence(seq, events, params);
    end
    
    % so seq and events now contain the starting point for the
    % optimization iterations
    
    % -------------- Run parameter-testing iterations ----------------
    for i=1:50
        
        % re-seed alignments, since seq changes sometimes
        events = seedaligns(seq,events,params);
        events = seedaligns(seq,events,params);
        events = seedaligns(seq,events,params);

        
        curscore = seqalign(refseq,seq);
        % display initial score
        if i==1
            fprintf('Max score: %0.2f\n',curscore);
        end
        
        fn = fieldnames(params);
        
        seqscores = zeros(nt,1);
        seqs = cell(nt,1);
        allparams = cell(nt,1);
        
        % now run parallel trials
        parfor j=1:nt
            curparams = params;
            
            % change some parameters by a few percent or so
            nc = randsubset(1:numel(fn),2);
            for k=nc
                curparams.(fn{k}) = curparams.(fn{k}) * (0.1*randn() + 1);
            end
            % except not these two
            curparams.stripe_width = 150;
            curparams.lik_offset = 4.5;
            
            % and now run the mutation
            newseq = MutateSequence(seq, events, curparams);
            seqscores(j) = seqalign(refseq,newseq);
            allparams{j} = curparams;
            seqs{j} = newseq;
        end
        
        % now save the best one from this iteration
        [maxscore,maxind] = max(seqscores);
        if maxscore > curscore
            seq = seqs{maxind};
            params = allparams{maxind};
        end
        fprintf('Max score: %0.2f\n',maxscore);
        % and save the best parameteres
        PoreParams.Save('CS_params.conf',params);
        
    end

end

