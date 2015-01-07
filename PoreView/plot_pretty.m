function plot_pretty(sigdata, trange, filter, channels)
%PLOT_PRETTY Makes a print-worthy plot
%   plot_pretty(sigdata, trange, filter, channels)
%   Takes a SignalData object and a time range as input.
%   Stephen Fleming, September 19, 2014

    % downsample to about 200000 points
    pt = 200000;
    n = round((trange(2)-trange(1))/sigdata.si); % # data points
    newSampleFreq = min([(pt/n)*(1/sigdata.si), filter*10, 1/sigdata.si]); % 5e6 points, or 10x oversampling, or full sampling, whichever is less
    for i=1:numel(channels)
        [data, ~] = downsample(sigdata,channels(i),trange,filter,newSampleFreq); % filter at 1kHz too
        d{i} = data;
        dmed{i} = medfilt1(data,51,1e5); % median filter
    end
    
    % Create a time axis
    time = linspace(trange(1),trange(2),numel(dmed{1}));
    
    % Get the reduced version of the unfiltered data as well
    raw = sigdata.getViewData(trange);
    
    % Plot!
    figure(2)
    clf(2)
    h = gcf;
    plot(raw(:,1),raw(:,2)*1000,'Color',[0.85,0.85,0.85])
    hold on
    if sigdata.nsigs>2 % we have pulses recorded
        [pulses,~,~] = pulse_analysis(sigdata,trange,[],[]); % get pulse timings
        line(repmat(pulses',1,2)',repmat(get(gca,'ylim'),numel(pulses),1)','Color',[0.9 1 0.9],'LineStyle','-') % vertical lines
    end
    color{1} = 'b';
    color{2} = 'r';
    for i = numel(channels):-1:1
        plot(time,dmed{i},'Color',color{i})
    end
    xlabel('Time (s)','FontSize',28)
    ylabel('Current (pA)','FontSize',28)
    name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];
    t = inputdlg('Enter title:','Input');
    title(t,'FontSize',24)
    annotation('textbox', [0.75 0.87 0 0], 'String', name, 'FontSize', 20);
    xlim([trange(1), trange(2)])
    ylim([min(0,min(dmed{1})-20),max(dmed{1})+20])
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
    set(h,'Position',[100 500 1000 400]) % size the figure
    set(gca,'FontSize',24)
    
end

