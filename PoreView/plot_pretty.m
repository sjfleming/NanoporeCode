function plot_pretty(sigdata, trange, filter, channels, t)
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
    
    % pulses
    if sigdata.nsigs>2
        display('Finding pulse timings...')
        % find possible pulses quickly using a derivative
        pds = 1000;
        [pulsedata, ~] = downsample(sigdata,4,trange,pds*10,pds);
        difference = diff(pulsedata);
        [~,candidates] = findpeaks(difference.*(difference>0),'MinPeakHeight',2,'MinPeakDist',50);
        candidates = trange(1) + candidates/pds; % convert from index to time
        % refine pulse timing
        pulses = [];
        for i = 1:numel(candidates)
            ind = (candidates(i)-10e-3)/sigdata.si;
            pulses(end+1) = sigdata.findNext(@(x) x(:,4)>1, ind);
            pulses(end) = pulses(end) * sigdata.si; % convert from index to time
        end
        display('Done.')
    end
    
    % Create a time axis
    time = linspace(trange(1),trange(2),numel(dmed{1}));
    %time = linspace(0,trange(2)-trange(1),numel(dmed{1}));
    
    % Get the reduced version of the unfiltered data as well
    raw = sigdata.getViewData(trange);
    
    % Plot!
    figure(2)
    clf(2)
    h = gcf;
    plot(raw(:,1),raw(:,2)*1000,'Color',[0.85,0.85,0.85])
    %plot(raw(:,1)-raw(1,1),raw(:,2)*1000,'Color',[0.85,0.85,0.85])
    hold on
    color{1} = 'b';
    color{2} = 'r';
    for i = numel(channels):-1:1
        plot(time,dmed{i},'Color',color{i})
    end
    xlabel('Time (s)','FontSize',28)
    ylabel('Current (pA)','FontSize',28)
    name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];
    
    % pulses
    if sigdata.nsigs>2 % we have pulses recorded
        line(repmat(pulses',1,2)',repmat(get(gca,'ylim'),numel(pulses),1)','Color',[0.8500    0.3250    0.0980],'LineStyle','-') % vertical lines
    end
    
    title(t,'FontSize',24)
    annotation('textbox', [0.75 0.87 0 0], 'String', name, 'FontSize', 20);
    %xlim([trange(1), trange(2)])
    xlim([0 trange(2)-trange(1)])
    ylim([min(0,min(dmed{1})-50),max(dmed{1})+50])
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
    set(h,'Position',[100 500 1000 400]) % size the figure
    set(gca,'FontSize',24)
    
%     figure(3)
%     clf(3)
%     h = gcf;
%     plot(time-time(1),dmed{1})
%     xlabel('Time (s)','FontSize',28)
%     ylabel('Current (pA)','FontSize',28)
%     xlim([0, trange(2)-trange(1)])
%     ylim([min(0,min(dmed{1})-10),max(dmed{1})+10])
%     set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
%     set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
%     set(h,'Position',[100 300 1000 400]) % size the figure
%     set(gca,'FontSize',24)
    
end

