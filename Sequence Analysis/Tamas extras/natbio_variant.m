function natbio_variant( trial )

    make all
    
    alparams = [];
    alparams.stripe_width = 100;
    alparams.insert_prob = 0.01;
    alparams.lik_offset = 4.5;
    alparams.do_fast = true;
    
    if nargin < 1
        trial = 1;
    end

    coverages = [1 2 4 8 15 25];
    ncoverages = numel(coverages);
    % -------------- Allocate output arrays ----------------
    variantscores = zeros(ncoverages,1);
    variantcounts = zeros(ncoverages,2);
    
    % -------------- Run de novo iterations ----------------
    for c=1:numel(coverages)
        % dunno why i'm doing this, random is random
        rseed = 1+12*trial;
        if trial>20
            % run all insertion trials on same dataset
            rseed = 1337;
        end
        
        % get the events
        % use actual sequence of M13
        % (seq here gets set to M13)
        [events,M13] = natbio_init(rseed,coverages(c));
        events = order_events(M13,events);
        % initialize the alignments using mapalign
        % (do it a few times, just for good measure)
        events = seedaligns(M13,events,alparams);
        events = seedaligns(M13,events,alparams);
        events = seedaligns(M13,events,alparams);
        
        % now do the actual variant scoring...
        % score all possible point mutations
        mutscores = ScoreMutations(M13,events,alparams);
        
        % figure out where events align
        nal = zeros(numel(M13),1);
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
        
        % and trim down mutscores
        mutscores = mutscores((i0*8):(i1*8));

        % remove duplicates, within a threshold
        % (these are assumed to be the same or thereabouts)
        thresh = 1e-8;
        mutuniq = unique(thresh*round(mutscores/thresh));

        fprintf('Point errors for trial %d, coverage %d: %0.2f (%0.2f, %d to %d)\n',trial,coverages(c),...
            100*sum(mutuniq<0)/numel(mutuniq),100*sum(mutscores<0)/numel(mutscores),i0,i1);
        
        variantscores(c) = 100*sum(mutuniq<0)/numel(mutuniq);
        variantcounts(c,:) = [sum(mutuniq<0),numel(mutuniq)];
        save(sprintf('./Out/Var_%03d.mat',trial),'variantscores','variantcounts','coverages');
    end        
end

