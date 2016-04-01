function plot_pretty(sigdata, trange, filter, channels, t)
%PLOT_PRETTY Makes a print-worthy plot
%   plot_pretty(sigdata, trange, filter, channels)
%   Takes a SignalData object and a time range as input.
%   Stephen Fleming, September 19, 2014

    if strcmp('none',filter)
        data = sigdata.getByTime(trange);
        for i = 1:numel(channels)
            if channels(i)==3 % voltage
                data(:,i) = data(:,channels(i))';
            else
                data(:,i) = data(:,channels(i))'*1000; % put in pA
            end
        end
    else
        % downsample to about 50000 points
        filtname = sprintf('Low-pass Bessel (%d Hz)', filter);
        f_chan = sigdata.addVirtualSignal(@(d) filt_lpb(d,4,filter), filtname);
        channels(1) = f_chan(1);
        pt = 5000;
        pt_actual = diff(trange)/sigdata.si; % number of points in raw data
        r = round(pt_actual/pt); % decimation factor
        data = [];
        approx = ceil(pt_actual/9e5);
        chnk = pt_actual/approx;
        for i=trange(1)/sigdata.si:chnk:(trange(2)/sigdata.si-chnk)
            d = [];
            d = sigdata.get(i:min(i+chnk,trange(2)/sigdata.si), channels);
            d(:,1) = d(:,1)*1000; % to pA
            %d = filt_lpb(d,4,filter);
            %data = [data; filt_decimate(d,r)];
            %dmed{i} = medfilt1(data,51,1e5); % median filter
            dfinal = [];
            if size(d,1)>2*pt/(pt_actual/chnk)
                for c = 1:numel(channels)
                    r = round(size(d,1)/(pt/(pt_actual/chnk)));
                    d1 = accumarray(1+floor((1:size(d,1))/r)',d(:,c)',[],@max);
                    d2 = accumarray(1+floor((1:size(d,1))/r)',d(:,c)',[],@min);
                    dfinal(:,c) = reshape([d1, d2]',1,[])';
                end
            end
            data = [data; dfinal];
        end
    end

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

    % Create a time axis
    time = linspace(trange(1),trange(2),size(data,1))';
    %time = linspace(0,trange(2)-trange(1),size(data,1))';

    % Get the reduced version of the unfiltered data as well
    raw = sigdata.getViewData(trange);
    % if view data is too long, reduce it
    if size(raw,1)>2*pt
        r = round(size(raw,1)/pt);
        raw1 = accumarray(1+floor((1:numel(raw(:,2)))/r)',raw(:,2)',[],@max);
        raw2 = accumarray(1+floor((1:numel(raw(:,2)))/r)',raw(:,2)',[],@min);
        clear raw
        raw(:,2) = reshape([raw1, raw2]',1,[]);
        raw(:,1) = linspace(trange(1),trange(2),size(raw,1));
    end
    
%     % Add an offset manually
%     m = 0.0065; % slope of offset, pA / second
%     t_init = trange(1); % seconds, where offset correction starts
%     offset = m*(time-t_init)';
%     dmed{1} = dmed{1} + offset;
%     offset_raw = m*(raw(:,1)-t_init) * 1e-3;
%     raw(:,2) = raw(:,2) + offset_raw;

    % Plot!
    figure(2)
    clf(2)
    h = gcf;
    plot(raw(:,1),abs(raw(:,2))*1000,'Color',[0.85,0.85,0.85])
    %plot(raw(:,1)-raw(1,1),raw(:,2)*1000,'Color',[0.85,0.85,0.85])
    hold on
    color{1} = 'b';
    color{2} = 'r';
    for i = numel(channels):-1:1
        plot(time,abs(data(:,i)),'Color',color{i})
    end
    xlabel('Time (s)','FontSize',28)
    ylabel('Current (pA)','FontSize',28)
    name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];

    % pulses
    if sigdata.nsigs>2 % we have pulses recorded
        line(repmat(pulses',1,2)',repmat(get(gca,'ylim'),numel(pulses),1)','Color',[0.8500    0.3250    0.0980],'LineStyle','-') % vertical lines
    end

    title(t,'FontSize',24)
    annotation('textbox', [0.75 0.87 0 0], 'String', [name ' ' num2str(filter) 'Hz'], 'FontSize', 20);
    xlim([trange(1), trange(2)])
    %xlim([0, trange(2)-trange(1)])
    ylim([min(0,min(abs(data(:,1)))-20),max(abs(data(:,1)))+20])
    set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
    set(gca,'OuterPosition',[0.01 0.01 0.98 0.98]) % fit everything in there
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

