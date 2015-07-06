function seq = testmut( )

    makeall    

    %pathname = 'C:\Minion\Lambda-73'; 
    pathname = 'C:\Minion\M13a';
    
    filename = [pathname '_2d.fasta.bam'];
    %filename = [pathname '_template.fasta.bam'];

    % first, load bam stuff into memory
    try
        h=baminfo(filename,'ScanDictionary',true);
    catch
        fprintf(2,'Not found: %s\n',filename);
    end
    bm = BioMap(filename, 'SelectReference', h.ScannedDictionary{1});

    %lambda = fastaread('.\References\Lambda_NEB.fasta');
    lambda = fastaread('.\References\M13mp18_EcoRI.fasta');
    lambda = lambda.Sequence;
    
    % ---------- Find aligned strand indices ---------------
    pd = PoreData(pathname);
    % get headers that align with the file
    refstart = 3000;
    refend = 4000;
    bminds = bm.getIndex(refstart,refend,'overlap',refend-refstart);
    headers = bm.Header(bminds);
    % and get their event indices
    %evinds = zeros(size(headers));
    evinds = pd.getIndex(headers);
    %for i=1:numel(headers)
    %    evinds(i) = pd.getIndex(headers{i});
    %end
    
    fprintf('Starting with %d strands...\n',numel(headers));
    
    
    % ---------- Find a seed sequence ---------------
    lstr = lambda(refstart:refend);
    % pick one of the strands to start with
    seq = [];
    alignscore = 0;
    curind = numel(evinds);
    % loop until we find one good enough
    while numel(seq) < 800 || alignscore < 75
        seq = pd.getSequence(evinds(curind),'2d');
        if isempty(seq)
            curind = curind - 1;
            continue
        end
        if (bm.Flag(bminds(curind)) > 0)
            seq = seqrcomplement(seq);
        end
        [~,p,~] = swalign(seq,lstr,'Alphabet','NT');
        % use only aligned bases
        seq = p(1,:);
        seq = seq(seq ~= '-');
        % score and start
        alignscore = seqalign(seq,lstr);
        curind = curind - 1;
    end
    fprintf('Starting with %d bases, initial score %0.1f\n',numel(seq),alignscore);
    
    % start with lambda and go down?
    %seq = lstr;
    

    % ---------- Read in strands and sequences ---------------
    
    
    % get the events, with proper ones reversed
    evinds(bm.Flag(bminds) > 0) = -evinds(bm.Flag(bminds)>0);
    evinds = evinds(randperm(numel(evinds)));
    events = pd.getEvents(evinds);

    % and the sequences
    seqs = pd.getSequences(evinds,'2d');

    
    for i=1:numel(seqs)
        % align each one with the current sequence
        [~,p] = swalign(seqs{i},lstr,'Alphabet','NT');
        % and keep only the overlapping region to speed up the mutation
        s = p(1,:);
        s = s(s ~= '-');
        seqs{i} = s;
    end
    
    % patchy code to make deviations not matter
    tic
    
    rng(13);
    
    % start with a "random" permutation of the events
    %allevents = events(randperm(numel(events)));
    allevents = events;
    
    % and only some of them
    events = allevents(1:50);
    
    % the current record
    alignrecord = 98.9;
    
    ser = display_score();
    ser.print(alignscore);
    
    alparams = [];
    alparams.stripe_width = 100;
    alparams.insert_prob = 0.0075;
    alparams.stay_prob = 0.01;
    alparams.extend_prob = 0.03;
    alparams.do_fast = true;
    
    fprintf('Doing full alignment...\n');
    [~,~,events] = align_likes(seq, events, alparams);

    % start timing
    tic;
    
    % ---------- Initial loop using template seqs ---------------
    %ss = randsubset(seqs,numel(events)/2);
    ss = {events.sequence};
    ss = ss(1:2:end);
    for n=1:5
        [seq,events] = MutateFast(seq, ss, events, alparams);
        
        alignscore = seqalign(lstr,seq);
        alignscore = 0.1*round(10*alignscore);
        fprintf('1D mutations %d: %0.1f | %d\n',n,alignscore,numel(seq));
        ser.print(alignscore);
    end
    
    %save('mutseq.mat','seq');
    %load('mutseq.mat');
    %load('bestseq.mat')
    
    %scores = align_likes(events, seq);
    % keep only the best events
    %events = events(scores > 500);
    
    % seed alignment
    %fprintf('Doing full alignment...\n');
    %[~,~,events] = align_likes(events, seq, alparams);

    % launch a viewer for looking at the sequence as we mutate it
    pv = PoreView();
    pv.setView([0.5 1000.5]);
    zoop = onCleanup(@() delete(pv.fig));
    
    % ---------- Main loop & refining ---------------
    
    n = 1;
    while 1
        %{
        if numel(events) == 27
            alparams.insert_prob = 0.001;
        elseif numel(events) == 37
            alparams.insert_prob = 0.001;
        lse
            alparams.insert_prob = 0.001;
        end
        %}
        % plot the stats as we go
        stats = align_stats(seq, events, lstr, alparams);
        plot_stats(seq, stats, events, pv);
        %pause
        
        % generate mutations
        seqs = viterbi_mut(seq, events);
        % cheat and use correct sequence?
        %seqs{end+1} = lstr;        
        % and then use those mutations for great justice
        [newseq, events] = MutateFast(seq, seqs, events, alparams);
        
        
        if mod(n,4) == 1
            fprintf('Doing point mutations...\n');
            % do full test, no cutoffs
            newseq = MutateSequence(newseq, events, alparams);
        end
        
        seq = newseq;
        
        % and measure scores
        seqscore = sum(align_likes(seq, events, alparams));
        lscore = sum(align_likes(lstr, events, alparams));

        alignscore = seqalign(lstr,seq);
        alignscore = 0.1*round(10*alignscore);
        fprintf('%d|%d/%d: %0.1f | %0.1f/%0.1f | %d\n',n,numel(events),numel(allevents),...
            alignscore,seqscore,lscore,numel(seq));
        if alignscore > alignrecord + 1e-3
            alignrecord = alignscore;
            system(sprintf('pythonw.exe ./Python/sendmsg.py %0.1f%%!',alignrecord));
        end
		%save(sprintf('%dpct.mat',floor(alignscore)),'seq');
		
        ser.print(alignscore);
        
        if n>=5
            % we have convergence, add a strand
            if numel(events) < numel(allevents)
                % add a few events
                nadd = min(10, numel(allevents)-numel(events));
                events(end+1:end+nadd) = allevents(numel(events)+1:numel(events)+nadd);
                n = 0;
            else
                break
            end
        end
        n = n + 1;
    end
    
    [s, p] = swalign(lambda, seq, 'Alphabet', 'NT')
    
    toc
end

