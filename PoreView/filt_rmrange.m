function d = filt_rmrange(d, range)
    %FILT_RMRANGE Removes ranges of points with their averages
    %   Data is passed with colums of [time, sig1, sig2, ...], and
    %       returns the same array but with the data filtered. Range should
    %       be an array [t0 t1 val; t0 t1 val; ...] where the points
    %       t0...t1 are replaced with val. Blame Ryan for this function.

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
            d(pts,j) = range(i,j+1);
        end
    end
end