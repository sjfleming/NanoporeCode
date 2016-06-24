function d_red = reduce_vec( d, nhalve )
% Trims a vector in min-max mode n times (halves n times)
    
    if (nhalve < 1)
        d_red = d;
        return
    end
    
    if mod(numel(d),2^(nhalve+1)) ~= 0
        d = d(1:end-mod(numel(d),2^(nhalve+1)));
    end

    nd = 2^(nhalve+1);
    % turn into array with each column stuff and stuff
    dr = reshape(d,[nd numel(d)/nd]);
    % now each column is nd adjacent elements, alternate
    % min and max...
    %[2*numel(d)/nd 1] --> [2*numel(d)/nd/size(d,2) size(d,2)]
    d_red = reshape([min(dr); max(dr)],[2*numel(d)/(nd*size(d,2)), size(d,2)]);
    % seriously, just trust me guys...

end

