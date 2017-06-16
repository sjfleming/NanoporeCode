function plot_pretty(sigdata, trange, filter, channels, t)
%PLOT_PRETTY Makes a print-worthy plot
%   plot_pretty(sigdata, trange, filter, channels)
%   Takes a SignalData object and a time range as input.
%   Stephen Fleming, September 19, 2014
    
    % filter
    filtname = sprintf('Low-pass Bessel (%d Hz)', filter);
    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_lpb(d,4,filter),filtname), sigdata, 'UniformOutput', false);
    fsigs = fsigs{1};
    
    % downsample
    data(:,2) = util.downsample_minmax(sigdata, fsigs(1), trange, 5000)'*1000;
    data(:,1) = linspace(trange(1),trange(2),size(data,1))';
    raw(:,2:3) = util.downsample_minmax(sigdata, 'view', trange, 5000)';
    raw(:,1) = linspace(trange(1),trange(2),size(raw,1))';

    % pulses
    if sigdata.nsigs>2
        display('Finding pulse timings...')
        % find possible pulses quickly using a derivative
        pds = 1000;
        pulsedata = util.downsample_pointwise(sigdata, 4, trange, 1e5);
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

    % Plot!
    figure(2)
    clf(2)
    h = gcf;
    tstart = data(1,1);
    plot(raw(:,1)-tstart,abs(raw(:,2))*1000,'Color',[0.85,0.85,0.85])
    hold on
    plot(raw(:,1)-tstart,raw(:,3),'Color','r')
    plot(data(:,1)-tstart,abs(data(:,2)),'Color','b')
    xlabel('Time (s)','FontSize',28)
    ylabel('Current (pA)','FontSize',28)
    try
        name = [sigdata.filename(end-27:end-20) '\_' sigdata.filename(end-7:end-4)];
    catch ex
        name = [];
    end

    % pulses
    if sigdata.nsigs>2 % we have pulses recorded
        line(repmat(pulses',1,2)'-tstart,repmat(get(gca,'ylim'),numel(pulses),1)','Color',[0.8500    0.3250    0.0980],'LineStyle','-') % vertical lines
    end

    title(t,'FontSize',24)
    annotation('textbox', [0.75 0.87 0 0], 'String', [name ' ' num2str(filter) 'Hz'], 'FontSize', 20);
    xlim([trange(1), trange(2)]-tstart)
    %xlim([0, trange(2)-trange(1)])
    ylim([min(0,min(abs(data(:,2)))-20),max(abs(data(:,2)))+20])
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
    set(h,'Position',[100 500 1000 400]) % size the figure
    set(gca,'FontSize',24)

end

