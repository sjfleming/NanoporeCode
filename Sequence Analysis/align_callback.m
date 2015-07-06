function align_callback( cf, seqs, sd0, sd1, key )
% Callback function for CrampFit multi-level alignment

    if strcmp(key.Key, 'f')
        t = mean(cf.getView());
        i0 = find(seqs(:,5)<t,1,'last')
        t0 = seqs(i0:i0+1,1);
        t1 = seqs(i0:i0+1,2);
        t0
        t1
        
        [df0, freqs0] = plotnoise(sd0,t0,2);
        [df1, freqs1] = plotnoise(sd1,t1,2);
        
        
        ss(1) = sqrt(2*sum(0.5*(df0(1:end-1)+df0(2:end)).*diff(freqs0)));
        ss(2) = sqrt(2*sum(0.5*(df1(1:end-1)+df1(2:end)).*diff(freqs1)));
        ss(3) = sqrt(2*mean(df0(freqs0<1e4))*1.2e4);
        ss(4) = sqrt(2*mean(df1(freqs1<1e4))*1.2e4);
        
        
        d0 = sd0.getByTime(t0(1),t0(2));
        d1 = sd1.getByTime(t1(1),t1(2));
        d0 = d0(:,2);
        d1 = d1(:,2);
        disp([numel(d0) numel(d1)])
        mm = prctile([d0;d1],[1 99]);
        
        ss(5) = std(d0);
        ss(6) = std(d1)
        
        bins = linspace(mm(1),mm(2),50);
        h0 = histc(d0,bins);
        h1 = histc(d1,bins);
        h0 = h0 / max(h0);
        h1 = h1 / max(h1);
        
        % do gaussian mixture fit!
        mus = zeros(2);
        pees = [0.5 0.5; 0.5 0.5];
        sig = min(ss);
        
        d = {d0, d1};
        
        for j=1:2
            mus(j,:) = prctile(d{j},[33 67]);
            for i=1:50
                % estimate weights
                wA = pees(j,1)*normpdf(d{j},mus(j,1),sig);
                wB = pees(j,2)*normpdf(d{j},mus(j,2),sig);
                % normalize
                wN = wA+wB;
                wA = wA./wN;
                wB = wB./wN;
                % recalculate
                pees(j,:) = [mean(wA) mean(wB)];
                mus(j,:) = [sum(wA.*d{j})/sum(wA), sum(wB.*d{j})/sum(wB)];
            end
        end
        % end mixture fit
        
        figure(3);
        %plot(bins,[h0,h1]);
        bins = linspace(mm(1),mm(2),100);
        d0 = sort(d0);
        d1 = sort(d1);
        f0 = cumsum(d0);
        f1 = cumsum(d1);
        [~, ia0] = unique(d0);
        [~, ia1] = unique(d1);
        d0 = d0(ia0);
        f0 = f0(ia0) / f0(end);
        d1 = d1(ia1);
        f1 = f1(ia1) / f1(end);
        f0 = interp1(d0,f0,bins);
        f1 = interp1(d1,f1,bins);
        
        f0 = conv(f0,0.2*[1 1 1 1 1],'same');
        f1 = conv(f1,0.2*[1 1 1 1 1],'same');
        f0([1:3, end-2:end]) = nan;
        f1([1:3, end-2:end]) = nan;
        
        bb = 0.5*(bins(2:end)+bins(1:end-1));
        db = bins(2)-bins(1);
        
        plot(bb,diff(f0)/db,bb,diff(f1)/db);
        hold on
        clrs = {'red','green'};
        for i=1:2
            for j=1:2
                plot(bb,pees(i,j)*normpdf(bb,mus(i,j),sig),clrs{i});
            end
            plot(bb,pees(i,1)*normpdf(bb,mus(i,1),sig)+pees(i,2)*normpdf(bb,mus(i,2),sig),clrs{i});
        end
    end
end

