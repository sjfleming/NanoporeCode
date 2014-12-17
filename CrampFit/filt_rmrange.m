function d = filt_rmrange(d, range)
    %FILT_RMRANGE Removes ranges of points
    %   Data is passed with colums of [time, sig1, sig2, ...], and
    %       returns the same array but with the data filtered.

    
    % this function should be given ranges as [t0 t1 val; t0 t1 val; ...]
    % (or [t0 t1 val0 val1] for files with multiple signals)

    % check every row of range
    for i=1:size(range,1)
        % do we have any points in the range?
        pts = and(d(:,1)>range(i,1),d(:,1)<range(i,2));
        if ~any(pts)
            continue
        end
        % loop through each column, and set all the points to values
        % in the last bits of the range array
        for j=2:size(d,2)
            d(pts,j) = NaN;%range(i,j+1);
        end
    end
end