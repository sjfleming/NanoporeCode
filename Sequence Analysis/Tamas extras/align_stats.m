function [stats,staystats] = align_stats(seq, events, refseq, params)
    % aight collect some ref_align stats and shit
    stats = zeros(numel(seq),8);
    staystats = zeros(numel(events),1);
    
    % alignment indices with reference/true sequence
    [~,p,st] = swalign(seq,refseq,'Alphabet','NT');
    % find gaps in our sequence, which are "errors"
    ginds = find(p(1,:) == '-');
    ginds = ginds(ginds>1)-1;
    % and make them into errors
    p(2,ginds) = ' ';
    % now remove all insertions to our sequence from p-column
    p = p(:,p(1,:) ~= '-');
    ginds = find(p(2,:) == '|') + st(1) - 1;
    ginds = ginds(ginds > 0);
    % and save the match locations n shit
    stats(ginds,1) = 1;
    
    % realign, and get ref likelihoods
    [~,events,seqreflike] = align_likes(seq,events,params);
    % get differential likelihoods
    stats(2:end,2) = diff(seqreflike);
    
    % and now do an accumulation loop
    for i=1:numel(events)
        % how many points at each location?
        ra0 = events(i).ref_align > 0;
        refcounts = accumarray(events(i).ref_align(ra0)+2,1,[numel(seq),1],@sum);
        stayn = accumarray(refcounts(refcounts>0),1,[],@sum);
        if numel(stayn) > numel(staystats)
            staystats(numel(staystats)+1:numel(stayn)) = 0;
        end
        staystats(1:numel(stayn)) = staystats(1:numel(stayn)) + stayn;
        % how many strands align?
        stats(:,3) = stats(:,3) + (refcounts > 0);
        % how many extra points/stays?
        stats(:,4) = stats(:,4) + (refcounts>1).*(refcounts-1);
        % contigs: have a point at ref-1 and ref and ref+1
        contigs = refcounts(3:end)>0 & refcounts(2:end-1)>0 & refcounts(1:end-2)>0;
        stats(2:end-1,5) = stats(2:end-1,5) + contigs;
        % check for insertions and put them on the prev state
        insinds = find(events(i).ref_align == -1) - 1;
        insinds = events(i).ref_align(insinds(insinds > 0));
        insinds = insinds(insinds > 0)+2;
        stats(insinds,6) = stats(insinds,6) + 1;
        % how many strands overall
        rrng = min(events(i).ref_align(ra0))+2:max(events(i).ref_align(ra0))+2;
        stats(rrng,7) = stats(rrng,7) + 1;
        % add up the times everywheres
        reftimes = accumarray(events(i).ref_align(ra0)+2,events(i).length(ra0),[numel(seq),1],@sum);
        %reftimes(reftimes > 0.1) = 0.1;
        stats(:,8) = stats(:,8) + reftimes;
        % and how many insertions came before this guy
    end
    
    %stats(:,8) = smooth_gauss(stats(:,8),2);
    stats(:,8) = medfilt1(stats(:,8),5);
    
    % trim staystats
    staystats = staystats(1:find(staystats==0,1,'first'));
end