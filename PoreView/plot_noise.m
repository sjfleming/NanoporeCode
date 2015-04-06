function plot_noise(sigdata, trange)
%PLOT_NOISE Makes or adds to existing noise plot
%   plot_noise(sigdata, trange)
%   Takes a SignalData object and a time range as input, uses the same
%       algorithm as ClampFit (I think).

    % do same thing as ClampFit does - average spectral segs
    % ClampFit does 2*65536 pts per seg
    % get start and end index
    irange = floor(trange/sigdata.si);
    
    fftsize = 2*65536;

    % only process real signals
    dfft = zeros(fftsize,sigdata.nsigs);
    % number of frames
    nframes = 0;
    
    wh = waitbar(0,'Calculating power spectrum...','Name','PoreView');
    
    for ind=irange(1):fftsize:irange(2)
        % get only the real signals
        d = sigdata.get(ind:ind+fftsize-1,1+(1:sigdata.nsigs));
        % quit if we don't have enough points
        if size(d,1) < fftsize
            break
        end
        % calculate power spectrum
        df = sigdata.si*abs(fft(d)).^2/fftsize;
        % and add to fft accum.
        dfft = dfft + df;
        nframes = nframes + 1;
        waitbar((ind-irange(1))/(irange(2)-irange(1)));
    end
    
    close(wh);

    % do the averaging
    dfft = dfft / nframes;
    
    % and calculate the frequency range
    freqs = 1/sigdata.si*(0:fftsize-1)/fftsize;
    
    % range of freqs to plot, only do half (Nyquist and all)
    imax = floor(size(dfft,1)/2);

    hf = findobj('Name','Noise Power Spectrum');
    if isempty(hf)
        % make a new plot figure, leave menu bar
        hf = figure('Name','Noise Power Spectrum','NumberTitle','off');

%         % Build a new color map
%         CO = [  0.0 0.2 0.7;...
%                 0.0 0.5 0.0;...
%                 0.7 0.1 0.1;...
%                 0.0 0.6 0.6;...
%                 0.5 0.0 0.5;...
%                 0.6 0.6 0.3  ];
% 
%         % Set the color order 
%         set(hf, 'DefaultAxesColorOrder', CO);
        
        % and make axes
        ax = axes('Parent',hf,'XScale','log','YScale','log',...
            'XLimMode','manual','YLimMode','auto','NextPlot','add',...
            'TickDir','out');
        
        % scale x axis
        set(ax,'XLim',[1 freqs(imax)]);
        set(gca,'FontSize',20)
        set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
        set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
        set(gcf,'Position',[100 500 800 400]) % size the figure
        % grid
        grid on
        grid minor
        box on
        % labels
        title('Noise Power Spectrum');
        ylabel('Power (nA^2/Hz)')
        xlabel('Frequency (Hz)')
    end

    % bring figure to front
    figure(hf);
    % and plot
    %plot(freqs(1:imax)',dfft(1:imax,:));
    plot(freqs(1:imax)',dfft(1:imax,1));
    name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];
    annotation('textbox', [0.75 0.9 0 0], 'String', name, 'FontSize', 20);
end

