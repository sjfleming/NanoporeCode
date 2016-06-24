function [ind1, ind2] = dtw_level_align(lev1, lev2)
% dtw_level_align(lev1,lev2) is a dynamic-time-warping level aligner
% lev1 is column vectors [level_means, level_stds] with level_stds optional
% same for lev2.
% dtw_level_align returns the indices of an aligned sequence, with ind1
% indexing the input lev1 and ind2 indexing lev2.
% there are no deletions in the alignment, only insertions.

% Stephen Fleming, 2016/06/21

    % cost function
    function x = cost(a, b)
        % this cost function determines which levels are aligned and how
        % much
        % a(1,1) is a level mean, and a(1,2) is an optional level std
        if numel(a)==1
            x = sqrt(abs(a(1,1)^2-b(1,1)^2));
        elseif numel(a)==2
            % give equal weight to the level mean and its noise (std)
            x = sqrt(abs(a(1,1)^2-b(1,1)^2)) + sqrt(abs(a(1,2)^2-b(1,2)^2));
        end
    end
    
    % initial array of indices
    DTW(1,1) = 0;
    DTW(2:size(lev1,1)+1,1) = Inf;
    DTW(1,2:size(lev2,1)+1) = Inf;
    
    % go through each array location and add up cost functions
    for i = 2:size(lev1,1)+1
        for j = 2:size(lev2,1)+1
            
            DTW(i,j) = cost(lev1(i-1,:),lev2(j-1,:)) + ...
                min([DTW(i-1,j), ...  % insertion
                DTW(i,j-1), ...       % deletion
                DTW(i-1,j-1)]);       % match
            
        end
    end
    
    % trace forward through the array to find the best aligned path
    % start at (2,2) effectively
    i = 1;
    j = 1;
    inds = [];
    while (i<size(DTW,1) && j<size(DTW,2))
        [~,next] = min([DTW(i+1,j), ...
            DTW(i+1,j+1), ...
            DTW(i,j+1)]);
        switch next
            case 1
                i = i+1;
            case 2
                i = i+1;
                j = j+1;
            case 3
                j = j+1;
        end
        inds(end+1,:) = [i,j];
    end
    
    % if we don't end at the end of both molecules, then put on some nans
    if inds(end,1)~=size(DTW,1)
        inds = [inds; [((1:(size(DTW,1)-inds(end,1))) + inds(end,1))', nan(size(DTW,1)-inds(end,1),1)]];
    elseif inds(end,2)~=size(DTW,2)
        inds = [inds; [nan(size(DTW,2)-inds(end,2),1), ((1:(size(DTW,2)-inds(end,2))) + inds(end,2))']];
    end
    
    % clean up so indices match levels passed in
    % not necessarily the same length if they don't align all the way...
    ind1 = inds(:,1)-1;
    ind1 = ind1(~isnan(ind1));
    ind2 = inds(:,2)-1;
    ind2 = ind2(~isnan(ind2));
    
end