function plot_errors(seq,events,mutscores,alparams)

    % make new figure
    fig = figure('Name','Event Info','NumberTitle','off');
    % and subplot axes on said figure
    hax1 = subplot(121);
    hax2 = subplot(122);
    
    % get the per-strand scores
    scores = align_likes(seq,events,alparams);
    
    % figure out order of mutation scores
    % first removing duplicates, within a threshold
    [mutuniq,ia] = unique(1e-5*round(1e5*mutscores));
    % ia is index of mutscores that are unique
    % we want to sort those, descending
    [~, inds] = sort(-mutuniq);
    errinds = ia(inds);    
    
    curind = 1;
    
    % re-plot error # ind
    function plotind(ind)
        curind = min(max(ind,1),size(errinds,1));
        % now make new sequence, with the error
        seqind = 1+floor((errinds(curind)-1)/8);
        bases = 'ACGT';
        bases = bases(bases ~= seq(seqind));
        % now check which one it is
        ind = mod(errinds(curind)-1,8);
        if ind == 0
            % deletion
            newseq = [seq(1:seqind-1) seq(seqind+1:end)];
        elseif ind < 4
            % mutation
            newseq = [seq(1:seqind-1) bases(ind) seq(seqind+1:end)];
        else
            bases = 'ACGT';
            newseq = [seq(1:seqind-1) bases(ind-3) seq(seqind:end)];
        end
        [newscores,newevents] = align_likes(newseq,events,alparams);
        % plot both panels
        plot_events(seq,events,seqind,hax1);
        plot_events(newseq,newevents,seqind,hax2);
        % and now plot all the score differences
        dscores = newscores - scores;
        xl1 = get(hax1,'XLim');
        for i=1:numel(events)
            s = 'T';
            if ~isempty(strfind(events(i).model.name,'compl'))
                s = 'C';
            end
            text(xl1(2),-i,sprintf('  %s: %+0.2f',s,dscores(i)),'HorizontalAlignment','left',...
                    'VerticalAlignment','middle','Parent',hax1,'FontSize',11);
        end
        % and re-set figure's key function to move left and right as well
        prevfun = get(fig,'KeyPressFcn');
        function keyfun(~,e)
            switch (e.Key)
                case 'leftarrow'
                    plotind(curind-1);
                case 'rightarrow'
                    plotind(curind+1);
                otherwise
                    % fall through to previous key function
                    prevfun([],e);
            end
        end
        set(fig,'KeyPressFcn',@keyfun);
        title(hax1,sprintf('Score #%d: %0.2f (%0.2f)',curind,sum(dscores),mutscores(errinds(curind))))
    end
    plotind(30);
end