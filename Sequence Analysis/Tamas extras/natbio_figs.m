function natbio_figs()

    load('Out_M13c/Seq_001.mat');
    data = nan(20,numel(coverages));

    for i=1:20
        load(sprintf('Out_M13c/Seq_%03d.mat',i));
        num = sum(sequencescores>0);
        data(i,1:num) = sequencescores(1:num);
    end
    
    % find min and max envelopes
    minscores = min(data);
    maxscores = max(data);
    meanscores = nanmean(data);
    
    xticks = log(coverages);
    yticks = 86:2:100;
    
    % and make the plot
    figure(1);
    clf
    subplot(211)
    
    %fill([1:numel(coverages) fliplr(1:numel(coverages))],yfun([minscores fliplr(maxscores)]),[1.0 0.8 0.75],'EdgeColor',0.5*[1 0.8 0.75]);
    hold on
    errorbar(xticks,meanscores,meanscores-minscores,maxscores-meanscores,'.-','Color',[0 0 0])
    fs = 12;
    xlabel('Coverage','FontSize',fs);
    ylabel('% Accuracy','FontSize',fs);
    set(gca,'FontSize',fs,'XTick',xticks,'XTickLabel',num2cell(coverages),'YTick',yticks,'YTickLabel',num2cell(yticks));
    dx = 0.3;
    xlim([xticks(1)-dx,xticks(end)+dx])
    ylim([min(yticks),max(yticks)]);
    box off
    for i=2:numel(yticks)-1
        h = line(xlim(),yticks(i)*[1 1]);
        set(h,'Color',0.7*[1 1 1],'LineStyle','-');
        uistack(h,'bottom');
    end
    pbaspect([2 1 1])
    
    % ----------------------------------------------------------
    
    load('Out_M13c/Var_001.mat');
    data = nan(20,numel(coverages));

    for i=1:20
        load(sprintf('Out_M13c/Var_%03d.mat',i));
        num = sum(variantscores>0);
        data(i,1:num) = variantscores(1:num);
    end
    
    % find min and max envelopes
    minscores = min(data);
    maxscores = max(data);
    meanscores = nanmean(data);
    
    % and yticks, and y-transform function
    xticks = log(coverages);
    yticks = 96.0:0.5:100;
    
    subplot(212)
    %fill([1:numel(coverages) fliplr(1:numel(coverages))],yfun([minscores fliplr(maxscores)]),[1.0 0.8 0.75],'EdgeColor',0.5*[1 0.8 0.75]);
    hold on
    errorbar(xticks,meanscores,meanscores-minscores,maxscores-meanscores,'.-','Color',[0 0 0])
    fs = 12;
    xlabel('Coverage','FontSize',fs);
    ylabel('% Accuracy','FontSize',fs);
    ylabels = num2cell(yticks);
    for i=2:2:numel(ylabels)
        ylabels{i} = '';
    end
    set(gca,'FontSize',fs,'XTick',xticks,'XTickLabel',num2cell(coverages),'YTick',yticks,'YTickLabel',ylabels);
    dx = 0.3;
    xlim([xticks(1)-dx,xticks(end)+dx])
    ylim([min(yticks),max(yticks)]);
    box off
    for i=2:numel(yticks)-1
        h = line(xlim(),yticks(i)*[1 1]);
        set(h,'Color',0.7*[1 1 1],'LineStyle','-');
        uistack(h,'bottom');
    end
    pbaspect([2 1 1])
    
    set(gcf,'Position',[100 100 600 700]);
    
    %% ----------------------------------------------------------
    
    % and the accuracy bar plot
    load lambda73
    figure(2)
    clf
    locs = {'Template','Complement','2D'};
    colors = [194 224 250; 210 226 139; 255 154 200]/255.0;
    xs = linspace(0,100,21);
    ys = 0*[xs; xs];
    for j=1:numel(locs)
        if isempty(idents{j})
            continue
        end
        
        ys(j,:) = hist(idents{j}(:,1), xs);
    end
    bs = bar(xs,ys','stacked');
    for j=1:numel(locs)
        set(bs(j),'FaceColor',colors(j,:));
    end
    ylabel('# of aligned strands','FontSize',13);
    xlabel('% identity','FontSize',13);
    set(gca,'FontSize',13);
    h=legend(locs);
    set(h,'FontSize',9);
    xlim([0 100])
    pbaspect([2 1 1])
    set(gcf,'Position',[500 200 750 400]);

    %% ----------------------------------------------------------
    
    figure(3);
    numx = 1000;
    numst = size(idents{3},1);
    % scale indices down to ref
    idents{3}(:,2:3) = numx*idents{3}(:,2:3)/max(idents{3}(:,3));
    pcov = zeros(numst,numx);
    for j=1:numst
        inds = round(idents{3}(j,2:3));
        inds = min(max(inds,1),numx);
        pcov(j,inds(1):inds(2)) = -1;
    end
    pcov = flipud(pcov);
    pcolor(pcov);
    shading flat
    colormap gray
    axis off
    box on
    caxis([-2.0 0]);
    pbaspect([1.5 1 1])
    set(gcf,'Position',[700 300 600 600]);
    
    
    %% ----------------------------------------------------------
    % now the states plot in fig1c
    
    figure(4);
    evs = natbio_init(2,1);
    ev = evs(1);

    [~,inds] = sort(ev.model.level_mean);
    lvls = ev.model.level_mean(inds);
    stdvs = ev.model.level_stdv(inds);
    kmers = ev.model.kmer(:,inds)';
    
    % draw the curve and filled part
    fill([(lvls-stdvs)' flipud(lvls+stdvs)'],[1:numel(inds) fliplr(1:numel(inds))],[1.0 0.8 0.75],'EdgeColor','none');%0.5*[1 0.8 0.75]);
    hold on
    plot(lvls,1:numel(lvls));
    xlabel('I (pA)');
    ylabel('5-mer')
    ylim([1 1024]);
    xlim([min(lvls) max(lvls)])
    title('')
    pbaspect([2 1 1])
    
    % now pick out the specific states to draw
    inds = [350 670 750 958 1010];
    
    for i=1:numel(inds)
        % first, draw a horizontal line
        plot([min(lvls) lvls(inds(i))],inds(i)*[1 1],'-');
        % and highlighted error bar line
        plot(lvls(inds(i))+stdvs(inds(i))*[-1 1],inds(i)*[1 1],'r-','LineWidth',1.5);
        % and label whee
        %text(min(lvls)+2.5,inds(i),kmers(inds(i),:),'VerticalAlignment','bottom','FontSize',8)
        text(lvls(inds(i))-4,inds(i),kmers(inds(i),:),'VerticalAlignment','top','FontSize',8)
        % and vertical line
        plot(lvls(inds(i))*[1 1],[inds(i) 1],'-');
        % and gaussian distribution
        mu = lvls(inds(i));
        sd = stdvs(inds(i));
        xs = linspace(mu-3*sd,mu+3*sd,200)';
        fill(xs,6e2*normpdf(xs,mu,sd),[194 224 250]/255)
        %alpha 0.5
    end
    set(gca,'YTick',[1 256 512 768 1024],'YTickLabel',num2cell([1 256 512 768 1024]));
    
    
    %% ----------------------------------------------------------
    
    % load a viewer and the data
    pv = pv_launch('C:\Minion\minion-pc_lambda_raw_5518_1.fast5');
    pd = PoreData('C:\Minion\Lambda-burnin');
    % find which indices match raw data
    inds = find(cellfun(@(x) any(strfind(x,'raw')>0),pd.Filename));
    % pick a channel
    ch = pd.Channel(inds(100));
    % and get all events in that channel
    pdinds = inds(pd.Channel(inds) == ch);
    % load it in PoreView
    pv.data = SignalData(pv.data.filename,'Channels',ch);
    pv.refresh();
    % and get the events
    events = pd.getEvents(pdinds);
    % and plot them
    hax = pv.getAxes(1);
    for i=1:numel(events)
        xs = doublemat(events(i).start);
        ys = doublemat(events(i).mean(1:end-1));
        plot(hax,xs(2:end-1),ys,'r','LineWidth',2);
    end
    
    %% ----------------------------------------------------------
    
    % generate accuracy figure or some shit
    pd = PoreData('C:\Minion\M13c');
    
    M13 = fastaread('./References/M13mp18_EcoRI.fasta');
    M13 = M13.Sequence;

    
    locs = {'t','c','2d'};
    for i=1:3
        evinds = find(pd.NumBases(:,i) > 6000 & pd.NumBases(:,i) < 8000);
        scores = 0*evinds;
        for j=1:numel(evinds)
            seq = pd.getSequence(evinds(j),locs{i});
            s1 = seqalign(M13,seq);
            scores(j) = max(s1,seqalign(M13,seqrcomplement(seq)));
            fprintf('%d/%d: %0.1f\n',j,numel(evinds),scores(j));
        end
        scatter(pd.NumBases(evinds,i),scores,'.')
    end
    
    %% ----------------------------------------------------------
    
    % and supplement figure with mutation finding ugh I'm getting tired of
    % this shit
    
    alparams = [];
    alparams.stripe_width = 150;
    alparams.insert_prob = 0.03;
    alparams.skip_prob = 0.04;
    alparams.stay_prob = 0.10;
    alparams.lik_offset = 4.5;
    alparams.do_fast = true;
    
    events = natbio_init(7,4);
    seq = events(3).sequence;
    events = order_events(seq,events);    
    events = seedaligns(seq,events,alparams);
    events = seedaligns(seq,events,alparams);
    events = seedaligns(seq,events,alparams);
    
    [~,~,reflike] = align_likes(seq,events,alparams);
    %
    seq1 = events(5).sequence;
    [~,seqal] = seqalign(seq,seq1);
    
    [~,~,seqlike] = align_likes(seq1,mapaligns(events,seqal),alparams);

    seqal = seqal-1;
    seqal = seqal(all(seqal>0,2),:);
    
    dlike = 0*seqal;
    dlike(:,1) = reflike(seqal(:,1));
    dlike(:,2) = seqlike(seqal(:,2));
    dlike = [0 0; diff(dlike)];
    dlikes = cusum(dlike(:,2) - dlike(:,1));
    
    figure;
    subplot(211)
    plot(diff([reflike(seqal(:,1)) seqlike(seqal(:,2))]))
    xlim([4055 4100])
    ylabel('Local Likelihood')
    subplot(212)
    plot(dlikes)
    xlim([4055 4100])
    xlabel('Sequence Position')
    ylabel('Cumulative Likelihood')


%% ----------------------------------------------------------

    % supplement figure on skip/stay params @ coverage 50
    
    events = natbio_init(1337,150);
    
    skips = arrayfun(@(x) x.model.skip_prob,events);
    stays = arrayfun(@(x) x.model.stay_prob,events);
    
    skips = reshape(skips,2,numel(skips)/2)';
    stays = reshape(stays,2,numel(stays)/2)';
    
    subplot(221)
    [n,c] = hist(skips);
    betterbar(c,n,0);
    xlim([0, 0.12])
    legend('Template','Complement')
    xlabel('Skip probability')
    ylabel('# of strands')
    
    
    subplot(222)
    [n,c] = hist(stays);
    betterbar(c,n,0);
    legend('Template','Complement')
    xlabel('Stay probability')
    
    
    % and now the best params data
    subplot(2,2,3:4)
        
    load('Out_M13c/Seq_001.mat');
    data = nan(40,numel(coverages));
    
    for i=1:40
        load(sprintf('Out_M13c/Seq_%03d.mat',i));
        num = sum(sequencescores>0);
        data(i,1:num) = sequencescores(1:num);
    end
    
    xticks = log(coverages);
    yticks = 89:2:99;
    
    hold on
    % find min and max envelopes
    %minscores = min(data(1:20,:));
    %maxscores = max(data(1:20,:));
    meanscores = nanmean(data(1:20,:));
    %errorbar(1:numel(coverages),yfun(meanscores),meanscores-minscores,maxscores-meanscores,'.-','Color',[1.0 0.5 0.45])
    plot(xticks,meanscores,'.-','Color',[1.0 0.5 0.45])

    minscores = min(data(21:40,:));
    maxscores = max(data(21:40,:));
    meanscores = nanmean(data(21:40,:));
    
    errorbar(xticks,meanscores,meanscores-minscores,maxscores-meanscores,'.-','Color',[0 0 0])
    
    legend('Estimated','Optimized','Location','SouthEast');
    
    fs = 12;
    xlabel('Coverage','FontSize',fs);
    ylabel('% Accuracy','FontSize',fs);
    set(gca,'FontSize',fs,'XTick',xticks,'XTickLabel',num2cell(coverages),'YTick',yticks,'YTickLabel',num2cell(yticks));
    
    dx = 0.3;
    xlim([xticks(1)-dx,xticks(end)+dx])
    ylim([88 100]);
    
    box off
    for i=1:numel(yticks)
        h = line(xlim(),yticks(i)*[1 1]);
        set(h,'Color',0.7*[1 1 1],'LineStyle','-');
        uistack(h,'bottom');
    end
    pbaspect([2 1 1])
    
end



