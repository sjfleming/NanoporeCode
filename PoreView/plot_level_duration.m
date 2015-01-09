function plot_level_duration(discreteData, name, tr, ti)
%PLOT_LEVEL_DURATION Makes a print-worthy histogram of level durations.
%   plot_level_duration(discreteData, V, name)
%   Takes a discreteData struct and the applied voltage and file name.

%   Stephen Fleming, September 24, 2014

    % make a plot of the level duration data
    figure(4)
    clf(4)
    h = gcf;
    t = (discreteData.levelTiming(:,2)-discreteData.levelTiming(:,1)); % level durations
    b = max(t)/75; % bin spacing in sec
    hist(t,b/2:b:(max(t)+b))
    xlabel('Duration (s)','FontSize',28)
    ylabel('Number of Levels','FontSize',28)
    title(['Level durations: ' ti],'FontSize',24)
    annotation('textbox', [0.72 0.85 0 0], 'String', [name  '\newline' num2str(tr(1),3) '-' num2str(tr(2),3) 's'], 'FontSize', 20);
    xlim([0 max(t)+2*b])
    %ylim([0,max(discreteData.levels)+10])
    set(h,'Position',[920 400 700 400]) % size the figure
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
    set(gca,'FontSize',24)
    
end