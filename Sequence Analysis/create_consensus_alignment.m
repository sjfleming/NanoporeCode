function consensus_levs = create_consensus_alignment(molecules)
% create_consensus_alignment(levs) does a best possible alignment of all
% the levels in levs using pairwise dynamic time warping and then keeping
% only the levels which show up in both pairs
% molecules is a cell array of molecule objects which need not be the same length

% Stephen Fleming 2016/06/21

    function consensus = pairwise_align(a,b)
        % align two sequences of levels pairwise
        
        % do the alignment using dynamic time warping
        [i1,i2] = dtw_level_align(a,b);
        
        % don't keep repeats in the alignment process
        [~,i1u,~] = unique(i1);
        [~,i2u,~] = unique(i2);
        num = max(numel(i1),numel(i2));
        levs = nan(num,2); % create space
        stds = nan(num,2);
        levs(i1u,1) = a(i1(i1u),1); % put in the aligned levels without repeats
        levs(i2u,2) = b(i2(i2u),1);
        stds(i1u,1) = a(i1(i1u),2);
        stds(i2u,2) = b(i2(i2u),2);
        consensus = [nanmean(levs,2), nanmean(stds,2)];
        
        % refine the consensus so that it tries to eliminate repeats
        event.level_means = consensus(:,1);
        event.level_stds = consensus(:,2);
        event.level_timing = ones(size(consensus));
        event.level_timing(:,2) = event.level_timing(:,1)+1;
        m = molecule(event);
        refined = m.get_robust_levels('pstay',-6,'pnoise',-10);
        clear consensus;
        consensus(:,1) = refined.level_means;
        consensus(:,2) = refined.level_stds;
        
        % find the points that agree really well
%         d = 20; % distance away in pA whicih is considered near-agreement
%         num = min(numel(i1),numel(i2));
%         levs = [a(i1(1:num),1),b(i2(1:num),1)];
%         stds = [a(i1(1:num),2),b(i2(1:num),2)];
%         consensus = [nanmean(levs(abs(levs(:,1)-levs(:,2))<d,:),2), nanmean(stds(abs(levs(:,1)-levs(:,2))<d,:),2)];
%         consensus = [mode(levs,2), mean(stds,2)];
        
    end

    % get levels from molecule objects
    levels = cell(1,numel(molecules));
    for i = 1:numel(molecules)
        % get levels reasonably free of noise and spurious stuff
        stdlim = 3.5; % levels with std greater than this are usually fake
        n = molecules{i}.get_robust_levels('pstay',-4,'pnoise',-10);
        levels{i} = [n.level_means(n.level_stds<stdlim), n.level_stds(n.level_stds<stdlim)];
    end
    
    % call pairwise_align and create consensus levels, then
    % keep doing this until you're down to one sequence of levels
    while numel(levels)>1
        % align random pairs
        pair = [true(1,2), false(1,numel(levels)-2)];
        pair = pair(randperm(numel(pair))); % logic for selecting the pair
        pair_levs = levels(pair);
        c = pairwise_align(pair_levs{1},pair_levs{2});
        % and replace those two sets of levels with the consensus
        levels = [levels(~pair), c];
    end
    
    consensus_levs = levels{1};

end