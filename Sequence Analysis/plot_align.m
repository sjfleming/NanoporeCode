function plot_align(events, refstates)
    % plots squiggle data
    
    pv = pv_launch();
    axes = pv.getAxes(1);
    
    % Build a new color map based on new Matlab colors
    CO = [       0    0.4470    0.7410;
            0.8500    0.3250    0.0980;
            0.9290    0.6940    0.1250;
            0.4940    0.1840    0.5560;
            0.4660    0.6740    0.1880;
            0.3010    0.7450    0.9330;
            0.6350    0.0780    0.1840];
    
    % and extend it
    CO = CO(1+mod(0:numel(events)-1,size(CO,1)),:);
        
    nref = numel(refstates);
    ts = (1:nref)';
    
    % reference squiggle plot
    p0 = plot(axes, ts, events(1).model.level_mean(refstates), '--', 'Color', [0 0 0], 'LineWidth', 2.0);
    
    % basic plot lines
    pp = zeros(numel(events),1);
    % full plot lines, showing skips and insertions
    pf = zeros(numel(events),1);
    
    function togglefocus(ind)
        % check if it's already got focus
        if get(pp(ind),'UserData')
            % set all to false, and to original colors
            for k=1:numel(pp)
                set(pp(k),'Color',CO(k,:),'UserData',0);
            end
            set(p0,'Color',[0 0 0]);
            set(pf, 'Visible', 'off');
        else
            % and make all the others light
            for k=1:numel(pp)
                set(pp(k),'Color',1-0.2*(1-CO(k,:)),'UserData',0);
            end
            % set it to true and dark
            set(pp(ind),'UserData',1,'Color',CO(1,:));
            set(p0,'Color',0.8*[1 1 1]);
            set(pf, 'Visible', 'off');
            set(pf(ind),'Visible','on');
            uistack(pp(ind),'top');
            uistack(pf(ind),'top');
        end
    end
    
    for i=1:numel(events)
        
        ra0 = events(i).ref_align > 0;
        % this accumulates all stays into a single point
        reflvls = accumarray(events(i).ref_align(ra0), events(i).mean(ra0), size(ts), @mean, nan);
        
        % now adjust for states
        reflvls = reflvls - events(i).model.level_mean(refstates) + events(1).model.level_mean(refstates);
        
        pp(i) = plot(axes, ts, reflvls, '.-', 'Color', CO(i,:), 'LineWidth', 1.0, 'UserData', 0);
        set(pp(i),'ButtonDownFcn',@(~,~) togglefocus(i));
        
        % and generate the more detailed plot, more slowly
        % first, create a time for each data point
        ts1 = zeros(numel(events(i).ref_align),1);
        lvls = ts1;
        curref = events(i).ref_align(1);
        curind = 1;
        while curind < numel(events(i).ref_align)
            % check if the next state is a real state
            if events(i).ref_align(curind) > 0
                curref = events(i).ref_align(curind);
                delta = -events(i).model.level_mean(refstates(curref)) + events(1).model.level_mean(refstates(curref));
                inds = find(events(i).ref_align == curref);
                if numel(inds) == 1
                    ts1(curind) = curref;
                    lvls(curind) = events(i).mean(curind) + delta;
                    curind = curind + 1;
                else
                    nr = numel(inds);
                    ts1(curind:curind+nr-1) = linspace(curref-0.1,curref+0.1,nr);
                    lvls(curind:curind+nr-1) = events(i).mean(curind:curind+nr-1) + delta;
                    curind = curind + nr;
                end
            else
                % if it isn't, assign point to just past last-used state,
                % and move on
                ts1(curind) = curref + 0.5;
                curind = curind + 1;
            end
        end
        pf(i) = scatter(ts1, lvls, 20, 'black', 'filled', 'Visible', 'off', 'HitTest', 'off');
    end
    
    pv.psigs(1).setY([min(rshape(reflvls)),max(rshape(reflvls))]);
    pv.setView([min(ts) max(ts)]);
    
end

