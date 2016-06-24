function analyze_errors()

    make all
    
    pd = PoreData('C:\Minion\M13c');
    
    M13 = fastaread('./References/M13mp18_EcoRI.fasta');
    M13 = M13.Sequence;
    
    % now we want to find all "useable" strands
    % these are defined as having a reasonable number of levels
    evinds = find(pd.NumBases(:,3) > 6000 & pd.NumBases(:,3) < 8000 & pd.NumEvents(:,2) > pd.NumEvents(:,1));
    
    fprintf('Found %d events total\n',numel(evinds));
        
    alparams = [];
    alparams.stripe_width = 150;
    alparams.insert_prob = 0.03;
    alparams.stay_prob = 0.10;
    alparams.extend_prob = 0.10;
    alparams.skip_prob = 0.05;
    alparams.lik_offset = 4.5;
    alparams.do_fast = true;

    
    events = pd.getEvents(evinds(1:20));

    seq = M13;

    % now find the orders of the events
    fprintf('Ordering eventssss');
    for j=1:numel(events)
        sfwd = seqalign(events(j).sequence(1:500),seq);
        srev = seqalign(seqrcomplement(events(j).sequence(1:500)),seq);
        fprintf('.');
        % is the strand probably backwards?
        if srev > sfwd
            events(j) = flip_event(events(j));
        end
    end
    fprintf('\n');
    
    events = clear_sdsd(events);
    % initialize the alignments using mapalign
    events = seedaligns(seq,events,alparams);
    events = seedaligns(seq,events,alparams);
    events = seedaligns(seq,events,alparams);

    % score all possible point mutations
    mutscores = ScoreMutations(seq,events,alparams);
    % and find unique ones
    mutuniq = unique(mutscores);
    fprintf('Point errors: %0.1f\n',100*sum(mutuniq<0)/numel(mutuniq));
    
    % and plot them both
    plot_errors(seq,events,mutscores,alparams);
end