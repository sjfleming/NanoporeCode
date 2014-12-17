function plot_signals(pv)
%PLOT_SIGNALS Makes a normal Matlab-style plot that copies currently
%   visible PoreView window, without all the UI junk
%   plot_signals(pv)
    
    % create the figure
    fig = figure('Name','PoreView Signals','NumberTitle','off');
    
    % set its position
    set(fig,'Units','normalized');
    set(fig,'Position',[0.1,0.1,0.8,0.8]);
    
    % set the axes color maps to be lighter if using reduced
    CO = get(pv.fig, 'DefaultAxesColorOrder');
    
    % now create subplots and copies of the panels and stuff
    npanels = numel(pv.psigs);
    
    hxax = pv.xaxes;
    
    if ~isempty(pv.data)
        signames = pv.data.getSignalList();
    end
    
    for i=1:npanels
        % copy the original axes from PoreView
        hax = copyobj(pv.getAxes(i),fig);
        % turn off grid lines, but keep ticks
        d = 0.015;
        set(hax,'Units','Normalized','OuterPosition',[d,(npanels-i)/npanels+d,1-2*d,1/npanels-2*d]);
        set(hax,'XTick',get(hxax,'XTick'),'YTick',get(pv.psigs(i).yaxes,'YTick'));
        set(hax,'XTickLabel',get(hxax,'XTickLabel'),'YTickLabel',get(pv.psigs(i).yaxes,'YTickLabel'));
        set(hax,'XGrid','off','YGrid','off','XColor',[0 0 0],'YColor',[0 0 0]);
        set(hax,'Box','on','XLimMode','manual','YLimMode','manual');
        set(hax,'LooseInset',[0 0 0 0]);
        
        
        % now manually plot some gridlines
        xt = get(hax,'XTick');
        yt = get(hax,'YTick');
        xl = get(hax,'XLim');
        yl = get(hax,'YLim');
        % Matlab is terrible
        xl = mean(xl) + 0.990*(xl-mean(xl));
        yl = mean(yl) + 0.950*(yl-mean(yl));
        
        for x=xt
            p = plot(hax,[x x],yl,'Color',0.92*[1 1 1],'Clipping','on');
            uistack(p,'bottom')
        end
        for y=yt
            p = plot(hax,xl,[y y],'Color',0.92*[1 1 1]);
            uistack(p,'bottom')
        end
        
        % set their labels and stuff
        if ~isempty(pv.data)
            title(hax,signames{pv.psigs(i).sigs(1)},'FontSize',14);
        else
            title(hax,'Plotted data','FontSize',14);
        end
        xlabel(hax,'Time (s)')
        ylabel(hax,'Current (nA)')
    end
    
    % also make a close callback, cause I'm spoiled and used to this
    function keyFun(~,e)
        if strcmp(e.Key,'escape')
            close(fig);
            return
        end
    end
    set(fig,'WindowKeyPressFcn',@keyFun);
    
end

