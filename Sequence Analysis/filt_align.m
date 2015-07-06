function [ d ] = filt_align( d, sd1, sd2, seqs )

    % total time in sequence 1 we're displaying
    ts1 = [d(1,1) d(end,1)];
    % indices of aligned sequences corresponding to times being displayed
    i0 = find(seqs(:,5)<ts1(1),1,'last');
    i1 = find(seqs(:,5)>ts1(2),1,'first');
    if (isempty(i0))
        i0 = 1;
    end
    if (isempty(i1))
        i1 = size(seqs,1);
    end
    % total time in sequence 2 we're displaying
    ts1 = [seqs(i0,1) seqs(i1,1)];
    ts2 = [seqs(i0,2) seqs(i1,2)];
    % data we're displaying, reduced version if too much
    if (ts1(2)-ts1(1)) > 6
        d1 = sd1.getViewData(ts1);
    else
        d1 = sd1.getByTime(ts1(1),ts1(2));
    end
    if (ts2(2)-ts2(1)) > 10
        d2 = sd2.getViewData(ts2);
    else
        d2 = sd2.getByTime(ts2(1),ts2(2));
    end
    % time-shift to query points, in both signals
    ts = interp1(seqs(i0:i1,5), seqs(i0:i1,1:2), d(:,1));
    % now set the signal
    d(:,2) = interp1(d1(:,1),d1(:,2),ts(:,1));
    d(:,3) = interp1(d2(:,1),d2(:,2),ts(:,2));
end

