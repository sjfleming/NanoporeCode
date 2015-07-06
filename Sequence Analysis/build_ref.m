function [ newref, dpath ] = build_ref(dpath, dref, d1)
%LEVEL_WARP LEVELS WTF AND STUDFF

    % if dpath doesn't contain every reference level, we want to augment it
    while dpath(1,1) > 1
        dpath = [dpath(1,1)-1, dpath(1,2); dpath];
    end
    % in both directions
    while dpath(end,1) < size(dref,1)+1
        dpath = [dpath; dpath(end,1)+1, dpath(end,2)];
    end

    % total number of levels
    k = size(dpath,1)-1;
    
    % reference sequence is array of size [#total levels, #of strands]
    % with NaNs where there are no levels
    
    newref = nan(k,size(dref,2)+1);

    % now set levels and times based on the destination index of the transition
    delts = dpath(2:end,:)-dpath(1:end-1,:);
    for i=1:k
        if delts(i,1) > 0
            % reference sequence
            newref(i,1:end-1) = dref(dpath(i,1),:);
        end
        if delts(i,2) > 0
            % new sequence
            newref(i,end) = d1(dpath(i,2));
        end
    end
end

