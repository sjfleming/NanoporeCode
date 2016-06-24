function [ pdf_out ] = plot_pdf( data )
%PLOT_PDF Plots sample pdf

    data = rshape(data);
    mm = prctile(data,[1 99]);
    bins = linspace(mm(1),mm(2),100);
    data = sort(data);
    f = cumsum(data);
    [~, ia] = unique(data);
    data = data(ia);
    f = f(ia) / f(end);
    f = interp1(data,f,bins);
    % smooth the pdf a bit
    f = conv(f,0.2*[1 1 1 1 1],'same');
    % and drop the ends
    f([1:3, end-2:end]) = nan;

    bb = 0.5*(bins(2:end)+bins(1:end-1));
    db = bins(2)-bins(1);

    plot(bb,diff(f)/db);
    
    pdf_out = [bb diff(f)/db];
end

