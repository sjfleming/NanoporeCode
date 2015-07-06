function [ seqs ] = time_warp( dpath, d0, d1, t0, t1)
%TIME_WARP WARP TIMES

    % add the ending time to the list to correspond to the time the 
    % last level ends, for the n+1th index

    k = size(dpath,1);
    ts = {t0,t1};
    ds = {d0,d1};

    % seqs has [t0 t1 lvl0 lvl1 tshift]
    % with NaNs where no level is drawn
    % each level line(s) at i should go from tshift(i) to tshift(i+1)
    seqs = nan(k,5);
    % bad times are zeros, not nans
    seqs(:,1:2) = 0;

    % now set levels and times based on the destination index of the transition
    delts = dpath(2:end,:)-dpath(1:end-1,:);
    for i=1:k-1
        for j=1:2
            seqs(i,j) = ts{j}(dpath(i,j));
            if delts(i,j) > 0
                seqs(i,j+2) = ds{j}(dpath(i,j));
            end
        end
    end
    seqs(end,1:2) = [t0(end) t1(end)];

    dts = sum(seqs(2:end,1:2) - seqs(1:end-1,1:2),2);
    dts = (t0(end)-t0(1))*dts/sum(dts);
    dts = [t0(1); dts];
    seqs(:,5) = cumsum(dts);

end

