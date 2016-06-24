function natbio_params()

    make all
    
    bestparams = [];
    if exist('bestparams.mat')
        load('bestparams.mat')
    else
        bestparams.insert = 0.02;
        bestparams.skip_t = 0.04;
        bestparams.skip_c = 0.04;
        bestparams.stay_t = 0.08;
        bestparams.stay_c = 0.08;
        bestparams.ext_t = 0.10;
        bestparams.ext_c = 0.10;
    end
    
    parc = parcluster('local');
    if isunix()
        parc.NumWorkers = 16;
    else
        parc.NumWorkers = 2;
    end
    nt = parc.NumWorkers;
    
    try
        parpool(parc);
    end
    
    % get the best sequence @ particular coverage as previously found
    load('Out/Seq_002.mat');
    curcov = 25;
    seq = sequences{find(coverages==curcov)};
    
    alparams = [];
    alparams.stripe_width = 100;
    alparams.insert_prob = 0.03;
    alparams.lik_offset = 4.5;
    alparams.do_fast = true;
    
    % -------------- Run parameter-testing iterations ----------------
    for i=1:10
        
        alparams.insert_prob = bestparams.insert;
        % get events from Seq_n
        [events,M13] = natbio_init(12*2 + 1,curcov);
        events = event_params(events);
        events = order_events(seq,events);
        % initialize the alignments using mapalign
        % (do it a few times, just for good measure)
        events = seedaligns(seq,events,alparams);
        events = seedaligns(seq,events,alparams);
        events = seedaligns(seq,events,alparams);
        
        %M13 = seqrcomplement(M13);
        curscore = seqalign(M13,seq);
        % display initial score
        if i==1
            fprintf('Max score: %0.2f\n',curscore);
        end
        
        fn = fieldnames(bestparams);
        
        % pick two params at random
        %{
        pinds = randi(numel(fn),[2 1]);
        f1 = fn{pinds(1)};
        f2 = fn{pinds(2)};
        % and generate stepped values for these params
        f1s = linspace(bestparams.(f1)*0.666,bestparams.(f1)*1.5,4);
        f2s = linspace(bestparams.(f2)*0.666,bestparams.(f2)*1.5,4);
        %}
        seqscores = zeros(nt,1);
        seqs = cell(nt,1);
        allparams = cell(nt,1);
        
        % now run parallel trials
        parfor j=1:nt
            curevents = events;
            curparams = bestparams;
            
            % change all parameters by +- 5% or so
            for k=1:numel(fn)
                curparams.(fn{k}) = curparams.(fn{k}) * (0.1*randn() + 1);
            end
            %curparams.(f1) = f1s(1+mod(j-1,4));
            %curparams.(f2) = f2s(1+floor((j-1)/4));
            alp = alparams;
            alp.insert_prob = curparams.insert;

            % set per-strand parameters
            for k=1:2:numel(events)
                % template strands
                curevents(k).model.skip_prob = curparams.skip_t;
                curevents(k).model.stay_prob = curparams.stay_t;
                curevents(k).model.extend_prob = curparams.ext_t;
                curevents(k+1).model.skip_prob = curparams.skip_c;
                curevents(k+1).model.stay_prob = curparams.stay_c;
                curevents(k+1).model.extend_prob = curparams.ext_c;
            end
            
            % and now run the mutation
            newseq = MutateSequence(seq, curevents, alp);
            seqscores(j) = seqalign(M13,newseq);
            allparams{j} = curparams;
            seqs{j} = newseq;
        end
        
        [maxscore,maxind] = max(seqscores);
        if maxscore > curscore
            seq = seqs{maxind};
            bestparams = allparams{maxind};
        end
        fprintf('Max score: %0.2f\n',maxscore);
        % and save the best parameteres
        save('bestparams.mat','bestparams');
        
    end

end

