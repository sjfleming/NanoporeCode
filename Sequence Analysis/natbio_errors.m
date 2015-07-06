function natbio_errors()

    % collect error stats from all output strands at a given coverage
    % (only really care about 50x)
    curcov = 50;
    [~,M13] = natbio_init();
    homotot = zeros(4,1);
    instot = zeros(4,1);
    deltot = zeros(4,1);
    muttot = zeros(4);
    for i=1:20
        load(sprintf('Out/Seq_%03d.mat',i));
        seq = sequences{coverages==curcov};
        if seqalign(seq,M13) < seqalign(seqrcomplement(seq),M13)
            seq = seqrcomplement(seq);
        end
        [homo,ins,del,mut] = error_stats(M13,seq);
        homotot = homotot + homo;
        instot = instot + ins;
        deltot = deltot + del;
        muttot = muttot + mut;
    end
    
    % now create the grid plot
    % rows are M13, need to transpose
    fullmat = [muttot deltot+homotot; instot' 0]';
    % change to %
    fullmat = fullmat*100/sum(sum(fullmat));

    figure();
    imagesc(fullmat,[-40 22])
    colormap gray
    pbaspect([1 1 1])
    
    for i=1:5
        for j=1:5
            text(j,i,sprintf('%0.1f',fullmat(i,j)),'HorizontalAlignment','center');
        end
    end
    
    set(gca,'XTick',1:5);
    set(gca,'YTick',1:5);
    set(gca,'XTickLabel',{'A','C','G','T','-'})
    set(gca,'YTickLabel',{'A','C','G','T','-'})
    set(gca,'TickLength',[0 0]);
    set(gca,'XAxisLocation','top')
    set(gca,'FontSize',16)
    ylabel('de novo')
    xlabel('M13mp18')

    % find homopolymer prevalence in M13
    bases = 'ACGT';
    hpcounts = zeros(4,1);
    for i=1:4
        hpinds = strfind(M13,repmat(bases(i),[1,5]));
        hpstr = M13;
        hpstr(:) = '0';
        for j=1:numel(hpinds)
            hpstr(hpinds(j):hpinds(j)+4)='1';
        end
        hpcounts(i) = sum(hpstr=='1');
    end
    
    colors = [0 114 189; 119 172 48; 77 190 238]/255.0;

    dtot = deltot + homotot;
    
    figure()
    subplot(211)
    bs = bar(100*[homotot deltot]/sum(dtot),0.65,'stacked');
    legend('Homopolymer','Non-homopolymer','Location','NorthWest')
    set(gca,'XtickLabel',{'A','C','G','T'})
    set(bs(1),'FaceColor',colors(2,:));
    set(bs(2),'FaceColor',colors(3,:));
    ylabel('% of all deletion errors')
    
    subplot(212)
    bs = bar(100*hpcounts/sum(hpcounts),0.65);
    set(gca,'XtickLabel',{'A','C','G','T'})
    ylabel('% of all homopolymer regions')
    set(bs,'FaceColor',colors(1,:));
end