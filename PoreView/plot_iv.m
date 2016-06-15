function [Vs, Is] = plot_iv(filename)
%PLOT_IV Plots IV curves from the file in filename, if it can.
%   IMPORTANT NOTE: Because I couldn't find the voltage information in the
%   abf file header, the voltages are hardcoded to be -200mV to 200mV, you
%   will need to modify if necessary.
%   Also, it only averages the last 25% of each sweep to get point, change
%   if you want (this is for very capacitive pores).

    try
        % try to load the entire file
        if filename(end-2) == 'a'
            [d,~,h]=abfload(filename);
        else
            [d,h]=cbfload(filename);
            % permute it to match abf
            d = permute(d,[2 3 1]);
        end
    catch
        fprintf(2,'Failed to load file %s as I-V!\n',filename);
        return
    end

    % check if it's an IV curve
    if (numel(size(d)) ~= 3)
        fprintf(2,'%s is not an IV curve.\n',filename);
        return
    end
    
    sz = size(d);
    
    % how much of the data do we want to use?
    % start with last 1/4, for now
    % indst = floor(sz(1)*0.75);
    
    % Get current values into a 1D/2D array, from 3D
    %Is = reshape(mean(d(indst:end,:,:),1),sz(2:3))';
    
    % use the mean of the second quarter after throwing out low current
    % blockages
    inds = 3000:121900;
    % Get current values into a 1D/2D array, from 3D
    Is = zeros(sz(3),1);
    d = reshape(d,sz(1),sz(3));
%     cut = max(abs(d(inds,:)),[],1)*96/100;
%     for i = 1:sz(3)
%         Is(i) = mean(d(abs(d(inds,i))>cut(i),i));
%     end
    
    figure(9)
    plot(d);
    nind = input('Input start and end points as a vector: ');
    for i = 1:sz(3)
        Is(i) = mean(d(nind(1):nind(2),i));
    end
    close(9);
    
    % make the voltages
    if Is(1)<0
        %Vs = linspace(-200,200,sz(3))';
        Vs = (-200:5:200)';
        p1 = 41;
        p2 = 51;
    else
        %Vs = linspace(200,-200,sz(3))';
        Vs = (200:-5:-200)';
        p1 = 31;
        p2 = 41;
    end
    if isfield(h,'setVoltages')
        Vs = h.setVoltages';
    end
    if numel(Is)<numel(Vs)
        display('File was truncated...')
        Is(end+1:numel(Vs)) = nan;
    end
    
    h = figure(8);
    clf(8)
    ym = max(abs(Is))+0.01;
    ylim([-1*ym ym])
    grid on
    box on
    line([0 0],[-1*ym ym],'Color','k')
    hold on
    line([-200 200],[0 0],'Color','k')
    plot(Vs,Is,'o:')
    
    hs = [];
    % fit the lines
    %coeffs = polyfit(Vs(1:ceil(numel(Is)/2)),Is(1:ceil(numel(Is)/2)),1);
    %p1 = floor(numel(Is)*1/4);
    %p2 = ceil(numel(Is)*2/4);
    coeffs = polyfit(Vs(p1:p2),Is(p1:p2),1);
    % draw the line
    ys = coeffs(1)*Vs + coeffs(2);
    fitLine = plot(Vs,ys,'LineStyle','--');
    legend(fitLine,sprintf('\\sigma = %4.4g nS',1000*coeffs(1)),'Location','northwest');
    name = [filename(65:68) '\_' filename(70:71) '\_' filename(73:74) '\_' filename(76:end-4)];
    annotation('textbox', [0.65 0.2 0 0], 'String', name, 'FontSize', 20);
    
    xlabel('Voltage (mV)','FontSize',28)
    ylabel('Current (nA)','FontSize',28)
    title('IV curve','FontSize',24)
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
    set(h,'Position',[100 500 600 500]) % size the figure
    set(gca,'FontSize',24)
    xlim([-200 200])
    
    assignin('base',['V_' filename(end-5:end-4)],Vs)
    assignin('base',['I_' filename(end-5:end-4)],Is)
    
%     hfig = figure('Name',sprintf('I-V Plot of %s',filename),'NumberTitle','off');
%     
%     % draw the axes
%     ax = axes('Position',[0.01 0.01 0.98 0.9],'Box','off','NextPlot','add');
%     % draw the data points
%     plts = plot(ax,Vs,Is,'Marker','o','MarkerSize',2,'MarkerEdgeColor','black');
%     % set them to filled circles
%     for i=1:length(plts)
%         set(plts(i),'MarkerFaceColor',get(plts(i),'Color'));
%     end
%     title('Quick I-V Plot');
%     % with lightly colored grid lines
%     set(ax,'XTickLabel','','YTickLabel','','XColor',0.8*[1 1 1],'YColor',0.8*[1 1 1],'FontSize',22);
%     grid on
%     % now manually make axes, cuz matlab is dumb
%     % first, extend axes limits slightly
%     set(ax,'XLimMode','manual','YLimMode','manual');
%     xlim = 1.1*get(ax,'XLim');
%     ylim = 1.1*get(ax,'YLim');
%     set(ax,'XLim',xlim);
%     set(ax,'YLim',ylim);
%     
%     % draw outer box
%     rectangle('Position',[xlim(1) ylim(1) (xlim(2)-xlim(1)), (ylim(2)-ylim(1))],...
%         'Clipping','off','EdgeColor','black');
% 
%     % and labels
%     dV = 40; % mV
% 	% V ticks, step from 20 to end-ish
%     Vts = dV/2:dV:Vs(end);
%     tickX = [-fliplr(Vts) Vts];
%     
%     % tick spacing in Y, do a clever thing
%     dy = ylim(2)/4;
%     ly = floor(log10(dy));
%     dy = dy/10^ly;
%     % dy is now between 1 and 10
%     % so pick spacing to be sensible numbers only (1, 2.5, 5, 10)
%     if (dy < 1.5)
%         dy = 1;
%     elseif (dy < 3.5)
%         dy = 2.5;
%     elseif (dy < 7.5)
%         dy = 5;
%     else
%         dy = 10;
%     end
%     dy = dy * 10^ly;
%     
%     tickY = dy:dy:ylim(2);
%     tickY = [-fliplr(tickY) tickY];
%     
%     % set the grid lines
%     set(ax,'XTick',tickX,'YTick',tickY);
%     
%     % now draw the labels
%     for i=1:length(tickX)
%         text(tickX(i),0,{'','',num2str(tickX(i))},'VerticalAlignment','middle',...
%             'HorizontalAlignment','center','FontSize',12,'Parent',ax);
%     end
%     for i=1:length(tickY)
%         text(0,tickY(i),['    ' num2str(tickY(i),'%#4.3g')],'VerticalAlignment','middle',...
%             'HorizontalAlignment','left','FontSize',12,'Parent',ax);
%     end
%     
%     % and add a bit to the lines
%     tickX = [2*tickX(1) tickX 2*tickX(end)];
%     tickY = [2*tickY(1) tickY 2*tickY(end)];
%     % draw axes ticks as black lines with '+' symbol
%     % super-duper hack wooooo
%     plot(ax,tickX,0*tickX,'black','Marker','+');
%     plot(ax,0*tickY,tickY,'black','Marker','+');
%     
%     % also display best-fit lines, while we're at it
%     xs = [tickX(1) tickX(end)];
%     hs = [];
%     for i=1:sz(2)
%         % fit the lines
%         coeffs = polyfit(Vs(floor(numel(Is)/2):end),Is(floor(numel(Is)/2):end,i),1);
%         % draw the line
%         ys = coeffs(1)*xs + coeffs(2);
%         hs(end+1) = plot(ax,xs,ys,'LineStyle','-.','Color',get(plts(i),'Color'));
%         % and set the legend
%         set(hs(end),'DisplayName',sprintf('\\sigma = %4.4g nS',1000*coeffs(1)));
%     end
%     
%     % now move original lines to the top
%     uistack(flipud(plts),'top');
%     
%     % get their names
%     M = get(hs,'DisplayName');
%     
%     % and display only those lines in the legend
%     legend(hs,M,'Location','NorthWest');
%     legend('show');
%     
%     % also make a close callback, cause I'm spoiled and used to this
%     function keyFun(~,e)
%         if strcmp(e.Key,'escape')
%             close(hfig);
%             return
%         end
%     end
%     set(hfig,'WindowKeyPressFcn',@keyFun);
end

