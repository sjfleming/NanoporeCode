function poisson_seq( trial )

    make all
    
    coverages = [2 4 8 15 25 50 100];
    ncoverages = numel(coverages);
    
    alparams = PoreParams.Load('CS_params.conf');
    
    % F5-proof
    if nargin < 1
        trial = 13;
    end
    
    [~,M13] = poisson_init(1,1);
    
    % calculates sequence score based on having maximal coverage
    function alignscore=calc_score(seq, events, typ, iter)
        % figure out where events align
        nal = zeros(numel(seq),1);
        for k=1:numel(events)
            ra0 = events(k).ref_align;
            ra0 = ra0(ra0 > 0);
            i0 = min(ra0);
            i1 = max(ra0);
            i0 = min(i0,numel(nal));
            i1 = min(i1,numel(nal));
            nal(i0:i1) = nal(i0:i1) + 1;
        end
        i0 = find(nal >= 0.75*numel(events),1,'first');
        i1 = find(nal >= 0.75*numel(events),1,'last');
        % now calculate best of forward and backward score
        subseq = seq(i0:i1);
        alignscore = seqalign(M13,subseq);
        alignscore = max(alignscore,seqalign(M13,seqrcomplement(subseq)));
        fullscore = seqalign(M13,seq);
        fullscore = max(fullscore,seqalign(M13,seqrcomplement(seq)));
        fprintf('%s mutations %d: %0.2f (%0.2f overall) | %d:%d / %d\n',typ,iter,alignscore,fullscore,i0,i1,numel(seq));
    end
    
    % -------------- Allocate output arrays ----------------
    sequencescores = zeros(ncoverages,1);
    sequences = cell(ncoverages,1);
    seqevents = cell(ncoverages,1);
    
    % -------------- Run de novo iterations ----------------
    for c=1:numel(coverages)
        
        outfile = sprintf('./Out/Seq_%03d.mat',trial);
        if exist(outfile)
            load(outfile);
            if sequencescores(c) > 0
                fprintf('Coverage %d found, skipping...\n',coverages(c));
                continue
            end
        end

        fprintf('de novo %d, coverage %d\n',trial,coverages(c));

        rseed = 1+12*trial;

        % get the events
        events = poisson_init(rseed,coverages(c));
        
        if c==1
            % pick the highest-scoring sequence from the starting (two) events
            bestlik = 0;
            bestseq = 0;
            for j=1:2:numel(events)
                seq = events(j).sequence;
                events = order_events(seq,events);
                events = seedaligns(seq,events,alparams);
                events = seedaligns(seq,events,alparams);
                liktot = sum(align_likes(seq,events,alparams));
                if liktot > bestlik
                    bestseq = seq;
                    bestlik = liktot;
                end
            end
            seq = bestseq;
        else
            % start with the previously-used sequence
            seq = sequences{c-1};
        end

        % order the events properly
        events = order_events(seq,events);

        % initialize the alignments using mapalign
        events = seedaligns(seq,events,alparams);

        % kk now mutate and optimize boom boom
        seqs = {};
        for n=1:numel(events)
            % don't save if it's already in there
            % which for 1/2 the strands, it will be
            if any(strcmp(seqs,events(n).sequence))
                continue
            end
            seqs{end+1} = events(n).sequence;
        end

        calc_score(seq,events,'Starting',0);

        n=1;
        while 1
            [seq,events,mutbases] = MutateFast(seq, randsubset(seqs,20), events, alparams);

            calc_score(seq,events,'1D',n);

            if mutbases < 5
                % done iterating
                break;
            end

            n = n + 1;
        end

        % and now start the Viterbi mutation trials
        % and do some point mutations, alternating
        for n=1:5
            seqs = viterbi_mut(seq, events);
            [seq,events,vmut] = MutateFast(seq, seqs, events, alparams);
            calc_score(seq,events,'Viterbi',n);

            [seq,events,pmut] = MutateSequence(seq, events, alparams);
            fprintf('Kept %d point mutations\n',pmut);
            calc_score(seq,events,'Point',n);

            if vmut+pmut < 10
                % done iterating
                break;
            end
        end

        finalscore = calc_score(seq,events,'Final',c);
        sequencescores(c) = finalscore;
        sequences{c} = seq;
        seqevents{c} = events;
        save(outfile,'sequencescores','sequences','coverages');
    end
    
end

