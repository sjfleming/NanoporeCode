function consensus_levs = create_consensus_alignment(levels)
% create_consensus_alignment(levels) does a best possible alignment of all
% the levels in levs using pairwise dynamic time warping and then keeping
% only the levels which show up in both pairs
% levels is a cell array, each of which is the matrix, in column vectors:
% [level_means, level_stds]

% Stephen Fleming 2017/08/03

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
        
        % find the points that agree really well
%         d = 0.001; % distance away in nA whicih is considered near-agreement
%         num = min(numel(i1),numel(i2));
%         levs = [a(i1(1:num),1),b(i2(1:num),1)];
%         stds = [a(i1(1:num),2),b(i2(1:num),2)];
%         consensus = [nanmean(levs(abs(levs(:,1)-levs(:,2))<d,:),2), nanmean(stds(abs(levs(:,1)-levs(:,2))<d,:),2)];
%         consensus = [mode(levs,2), mean(stds,2)];
        
        % eliminate repeats
        lev_diffs = abs(consensus(2:end,1) - consensus(1:end-1,1));
        std_compare = min(consensus(1:end-1,2), consensus(2:end,2));
        not_really_different = lev_diffs < std_compare; % neighboring levels are same within uncertainty
        not_really_different = [0; not_really_different];
        newcon = consensus(1,:);
        new_lev_inds = find(~not_really_different);
        for i = 1:numel(new_lev_inds)-1 % for all the ones that are new levels (i.e. not_really_different == 0)
            newcon(i+1,:) = [mean(consensus(new_lev_inds(i):new_lev_inds(i+1)-1,1)), mean(consensus(new_lev_inds(i):new_lev_inds(i+1)-1,2))];
        end
        newcon(i+1,:) = [mean(consensus(new_lev_inds(i):end,1)), mean(consensus(new_lev_inds(i):end,2))];
        consensus = cell(1);
        consensus{1} = newcon;
        
    end
    
    % call pairwise_align and create consensus levels, then
    % keep doing this until you're down to one sequence of levels
    while numel(levels)>1
        % align random pairs
        pair = [true(1,2), false(1,numel(levels)-2)]';
        pair = pair(randperm(numel(pair))); % logic for selecting the pair
        pair_levs = levels(pair);
        c = pairwise_align(pair_levs{1},pair_levs{2});
        % and replace those two sets of levels with the consensus
        inds = find(pair);
        levels = [levels(~pair); c];
    end
    
    consensus_levs = levels{1};

end