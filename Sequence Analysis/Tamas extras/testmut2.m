function str = testmut2( )
% test mutations on a single event/complement pair

    try
        mex align_like.cpp
        mex viterbikd.cpp
    end

    pathname = 'C:\Minion\Lambda-freezein';
    
    lambda = fastaread('.\References\Lambda_NEB.fasta');
    lambda = lambda.Sequence;
    
    pd = PoreData(pathname);

    evind = 131;
    
    % load these events
	[events(1), models(1)] = pd.getEvent(evind,'t');
	[events(2), models(2)] = pd.getEvent(evind,'c');
    % flip the complement side over
	events(2) = flip_event(events(2));
	models(2) = flip_model(models(2));
        
    % patchy code to make deviations not matter
    for i=1:0%numel(events)
        foo = 0*models(i).sd_stdv;
        models(i).sd_stdv = foo + 1/sqrt(2*pi);
        models(i).sd_mean = foo + 0;
        events(i).stdv = 0*events(i).stdv;
        % and let's make all the deviations the same while we're at it
        %models(i).level_stdv = foo + mean(models(i).level_stdv);
    end
    
    % run 1d viterbi to get template and complement sequences
    
    events(1).ref_align = (1:numel(events(1).ref_align))';
    events(2).ref_align = (1:numel(events(2).ref_align))';
    params = [];
    params.skip_prob = 0.1;
    params.stay_prob = 0.03;
    params.mutations = 4;
    params.mut_min = 0.5;
    params.mut_max = 1.0;
    
    seqs = {};
    seqs = [seqs statestoseq(viterbikd(events(1),models(1),params))];
    seqs = [seqs statestoseq(viterbikd(events(2),models(2),params))];
    
    alparams = [];
    alparams.stripe_width = 200;
    alparams.lik_offset = 5.0;
    
    % and get a few paths from both together...?
    [~,~,events] = align_likes(clearaligns(events),models,seqs{1},alparams);
    seqs = [seqs statestoseq(viterbikd(events,models,params))];
    [~,~,events] = align_likes(clearaligns(events),models,seqs{6},alparams);
    seqs = [seqs statestoseq(viterbikd(events,models,params))];
    
    % the sequence to iterate with
    seq = seqs{11};
    
    tic
    
    ser = display_score;
    
    pv = PoreView();
    pv.setView([0.5 7000.5]);
    zoop = onCleanup(@() delete(pv.fig));
    
    [~,p] = swalign(lambda,seqs{1},'Alphabet','NT');
    lstr = p(1,:);
    lstr = lstr(lstr ~= '-');

    events = clearaligns(events);
    
    % first run, just use existing mutations and don't generate new ones
    for n = 1:50
        stats = align_stats(seq, events, models, lstr, alparams);
        plot_stats(seq, stats, events, pv);
        
        [seq,events] = mutate_seq(seq, seqs, events, models, alparams);
        
        alignscore = seqalign(lambda,seq);
        ser.print(alignscore);
        fprintf('%d: %0.1f | %d\n',n,seqalign(lambda,seq),numel(seq));
        
        params.mutations = 10;
        [~,~,events] = align_likes(clearaligns(events),models,seq,alparams);
        seqs = statestoseq(viterbikd(events,models,params));
    end
    
end

