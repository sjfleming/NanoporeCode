function seqs = viterbi_mut(seq, events)

    % first, generate mutations using Viterbi
    vitev = randsubset(events,4);

    params = [];
    % first, generate some mutated paths
    params.mutations = 10;
    params.skip_prob = 0.05;
    dpaths = viterbikd(vitev,params);
    
    seqs = {};
    
    % then, use those to generate more paths! oh man such crazy wow
    for n=1:size(dpaths,2)
        seqs{end+1} = statestoseq(dpaths(:,n));
        continue;
        % get some new randomness up in here for each random thing
        vitind = randperm(numel(events),8);
        vitev = events(vitind);
        % and use that
        [~,~,newevents] = align_likes(vitev,seqs{end});
        newdpaths = viterbikd(newevents,params);
        % but skip the viterbi path one, cuz, like, who cares
        for i=2:size(newdpaths,2)
            seqs{end+1} = statestoseq(newdpaths(:,i));
        end
    end
end
