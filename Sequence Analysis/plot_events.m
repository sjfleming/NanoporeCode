function plot_events(seq,events,refind,hax)
% plot all of the events at a given reference index and stuff
    
    if nargin<4
        fig = figure('Name','Event Info','NumberTitle','off');
        clf(fig);
        hax = axes();
    else
        fig = get(hax,'Parent');
    end
    cla(hax);
    hold(hax,'on');
    
    % how many to pad to each side of refind
    npad = 8;
    seqind = [max(3,refind-npad),min(numel(seq)-2,refind+npad)];
    
    % cell array to hold indices corresponding to each reference position
    evinds = cell(numel(events),seqind(2)-seqind(1)+1);
    % and the insertions after each reference position
    insinds = cell(numel(events),seqind(2)-seqind(1)+1);
    
    % draw shaded rectangle for highlighted refind
    rectangle('Position',[refind-0.5,-numel(events)-0.5,1,numel(events)+1],...
        'Parent',hax,'FaceColor',1-0.1*[0 1 1],'EdgeColor','none');
    % now plot the sequence up top
    text(seqind(1):seqind(2),0.0*(seqind(1):seqind(2)),num2cell(seq(seqind(1):seqind(2))),...
        'HorizontalAlignment','Center','Parent',hax,'FontSize',14);
    % and vertical dotted lines
    for j=seqind(1):seqind(2)+1
        plot(hax,[j j]-0.5,[0.5 -numel(events)-0.5],':','Color',0.7*[1 1 1]);
    end
    % and horizontal lines
    for i=0:numel(events)
        plot(hax,[seqind(1)-0.5,seqind(2)+0.5],-[i i]-0.5,'-','Color',0.7*[1 1 1]);
    end
    
    set(hax,'XLim',[seqind(1)-0.5,seqind(2)+0.5]);
    
    % get the states, for later use
    states = seqtostates(seq);
    
    % function to move up and down when looking at events
    function keyfun(~,e)
        ind = get(fig,'UserData');
        haxs = get(fig,'Children');
        switch (e.Key)
            case 'uparrow'
                if isempty(ind)
                    ind = 2;
                end
                ind = ind - 1;
                set(haxs,'ylim',[-ind-1.5,-ind+1.5]);
            case 'downarrow'
                if isempty(ind)
                    ind = 0;
                end
                ind = ind + 1;
                set(haxs,'ylim',[-ind-1.5,-ind+1.5]);
        end
        set(fig,'UserData',ind);
    end
    set(fig,'KeyPressFcn',@keyfun);
    set(fig,'UserData',[]);
    zz = [];
    zz.Key='downarrow';
    keyfun([],zz);
    
    for i=1:numel(events)
        % populate cell arrays with indices
        for j=seqind(1):seqind(2)
            % event inds, just where ref_align is equal to that, except
            % minus 2 to account for state/sequence indexing
            evinds{i,j} = find(events(i).ref_align == (j-2));
            % insertion indices, where it's right after that state
            % and before the next one
            evafter = find(events(i).ref_align > (j-2),1,'first');
            if ~isempty(evafter)
                insinds{i,j} = find(events(i).ref_align(min(evinds{i,j}):evafter)==-1)+min(evinds{i,j})-1;
            end
        end
        mm = [min(events(i).model.level_mean)-5,max(events(i).model.level_mean)+5];
        % ok, now draw levels
        for j=seqind(1):seqind(2)
            % first, plot model state as a rectangle
            sd = events(i).model.level_stdv(states(j-2));
            ry = (events(i).model.level_mean(states(j-2))-mm(1)-sd)/diff(mm);
            rectangle('Position',[j-0.5,ry-i-0.5,1,2*sd/diff(mm)],...
                'Parent',hax,'FaceColor',0.92*[1 1 1],'EdgeColor','none');
            
            % if anything to draw
            if ~isempty(evinds{i,j})
                % find x indices
                xs = linspace(j-0.5,j+0.5,numel(evinds{i,j})+1)';
                % and y indices, scaled
                ys = (events(i).mean(evinds{i,j})-mm(1))/diff(mm);
                xs = doublemat(xs);
                ys = doublemat(ys)-i-0.5;
                plot(hax,xs(2:end-1),ys);
                text(mean(xs),-i-0.5,sprintf('%0.2f',diff(events(i).ref_like([evinds{i,j}(1)-1,evinds{i,j}(end)]))),...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','Parent',hax,'FontSize',8);
            end
            if ~isempty(insinds{i,j})
                % plot between levels
                ys = (events(i).mean(insinds{i,j})-mm(1))/diff(mm)-i-0.5;
                scatter(hax,0*ys + j + 0.5, ys, '.');
            end
        end
    end
    
end
