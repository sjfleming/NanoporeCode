function plot_squiggles(discreteData, name, tr, t)
%PLOT_SQUIGGLES Makes a print-worthy plot of level data at each base
%   plot_squiggles(discreteData, V, name)
%   Takes a discreteData struct and the applied voltage and file name.

%   Stephen Fleming, September 24, 2014

    % make a plot of the squiggle data
    figure(3)
    clf(3)
    h = gcf;
    plot(discreteData.levels,'o')
    hold on
    plot(discreteData.levels,'x')
    xlabel('Level','FontSize',28)
    ylabel('Current (pA)','FontSize',28)
    title(['Levels: ' t],'FontSize',24)
    annotation('textbox', [0.72 0.85 0 0], 'String', [name  '\newline' num2str(tr(1),3) '-' num2str(tr(2),3) 's'], 'FontSize', 20);
    xlim([0 numel(discreteData.levels)+1])
    ylim([0,max(discreteData.levels)+10])
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
    set(h,'Position',[200 400 700 400]) % size the figure
    set(gca,'FontSize',24)

end