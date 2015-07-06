function seq = laszlorun()

    makeall

    levents = [];
    % load the events file
    load laszloevents
    % load quadrometer map
    load laszloqmer
    % and the phix sequence
    phix = fastaread('.\References\phix174.fasta');
    phix = phix.Sequence;
    
    % convert the model into the oxford-type struct, by making all 4 5mers
    % corresponding to a 4mer identical
    pd = PoreData('C:\Minion\Lambda-freezein');
    event = pd.getEvent(1,'t');
    model = event.model;
    % kmer, variant, nothing happens
    model.weight = 1+0*model.weight;
    model.sd_mean = 0*model.sd_mean;
    model.sd_stdv = sqrt(0.5)+0*model.sd_stdv;
    qmers = [qmerdatabase.mean];
    qerrs = [qmerdatabase.error];
    inds = doublemat(doublemat((1:256)'));
    model.level_mean = 117*qmers(inds)';
    model.level_stdv = 250*qerrs(inds)';
    
    % now we need to format the events as events structs we're used to
    events = [];
    
    vparams = [];
    vparams.mutations = 0;
    vparams.skip_prob = 0.10;
    vparams.stay_prob = 0.05;
    
    alparams = [];
    alparams.stripe_width = 200;
    alparams.comp_cutoff = 0;
    alparams.skip_prob = 0.10;
    alparams.insert_prob = 0.10;
    alparams.stay_prob = 0.10;
    alparams.extend_prob = 0.20;
    alparams.offset = 5.5;

    
    pcts = [1 99];
    
    qr = prctile(model.level_mean,pcts);
    
    phix = seqreverse(phix);

    for i=1:numel(levents)
        if levents(i).count < 300
            continue
        end
        event.mean = rshape(levents(i).mean);
        e0 = event.mean > 0;
        event.mean = event.mean(e0);
        er = prctile(event.mean,pcts);
        event.mean = qr(1)+(event.mean-er(1))*diff(qr)/diff(er);
        %event.stdv = rshape(levents(i).stdv(e0));
        event.stdv = 0*event.mean;
        event.start = rshape(levents(i).start(e0));
        event.length = rshape(0*levents(i).mean);
        event.ref_align = rshape(1:numel(event.mean));
        event.ref_like = 0*event.ref_align;
        event.model = model;
        % and fuck it, run Viterbi on it to get each one some
        % sort of a sequence, y knot?
        event.sequence = statestoseq(viterbikd(event,vparams));
        event = clearaligns(event);
        [sfwd,evfwd] = align_like_omg(event,phix,alparams);
        [srev,evrev] = align_like_omg(event,seqrcomplement(phix),alparams);
        fprintf('Score: %f | %f\n',sfwd,srev);
        clf
        cdfplot(event.mean);hold on;cdfplot(event.model.level_mean)
        if isempty(events)
            events = event;
        else
            events(end+1) = event;
        end
    end
    
end

% function to merge nearby levels, for scaling
function levels=trimlevels(levels)
    %levels = sort(levels);
    while 1
        ds = abs(diff(levels));
        [m,ind] = min(ds);
        if m>1.5
            break;
        end
        levels(ind) = [];
    end
end