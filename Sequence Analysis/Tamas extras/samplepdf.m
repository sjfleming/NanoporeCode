function [ samppdf ] = samplepdf( data )
%SAMPLEPDF Get sample pdf

    mm = [min(data) max(data)];
    
    data = sort(rshape(data));

    % take cumulative sum
    f = cumsum(data);
    % drop repeated elements, from (sorted) data and cumulative sum
    [~, ia] = unique(data);
    data = data(ia);
    f = f(ia) / f(end);

    % and interpolate to our nice function
    bins = linspace(mm(1),mm(2),30);
    f = interp1(data,f,bins);
    % smooth the pdf a bit
    f = [0 0 f 1 1];
    f = conv(f,0.2*[0.5 1 2 1 0.5]','valid');
    % and drop the ends
    %f([1:3, end-2:end]) = nan;

    bb = 0.5*(bins(2:end)+bins(1:end-1));
    db = bins(2)-bins(1);

    samppdf = [bb; diff(f)/db];
    
end

