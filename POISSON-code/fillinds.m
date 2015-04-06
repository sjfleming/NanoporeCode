function inds = fillinds(inds)
% goes through and fills any zero or negative indices with previous values
% eg. [1 2 3 0 0 4 0 7 0 1 0] -> [1 2 3 3 3 4 4 7 7 1 1]

    i0 = inds(1,:);
    
    for i=2:size(inds,1)
        i0(inds(i,:) > 0) = inds(i,inds(i,:) > 0);
        inds(i,:) = i0;
    end

end