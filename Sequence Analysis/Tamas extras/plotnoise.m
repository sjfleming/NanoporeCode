function [dfft, freqs] = plotnoise(sigdata, trange, sigs)

    irange = round(trange/sigdata.si);    
    fftsize = irange(2)-irange(1)+1;

    d = sigdata.get(irange(1):irange(2),sigs);
    % calculate power spectrum
    dfft = sigdata.si*abs(fft(d)).^2/fftsize;
    
    % and calculate the frequency range
    freqs = 1/sigdata.si*(0:fftsize-1)/fftsize;
    
    % range of freqs to plot, only do half (Nyquist and all)
    imax = floor(size(dfft,1)/2);

    hf = findobj('Name','Noise Power Spectrum');
    if isempty(hf)
        % make a new plot figure, leave menu bar
        hf = figure('Name','Noise Power Spectrum','NumberTitle','off');

        % Build a new color map
        CO = [  0.0 0.2 0.7;...
                0.0 0.5 0.0;...
                0.7 0.1 0.1;...
                0.0 0.6 0.6;...
                0.5 0.0 0.5;...
                0.6 0.6 0.3  ];

        % Set the color order 
        set(hf, 'DefaultAxesColorOrder', CO);
        
        % and make axes
        ax = axes('Parent',hf,'XScale','log','YScale','log',...
            'XLimMode','manual','YLimMode','auto','NextPlot','add',...
            'TickDir','out');
        
        % scale x axis
        set(ax,'XLim',[1 freqs(imax)]);
        % grid
        grid on
        grid minor
        % labels
        title('Noise Power Spectrum');
        ylabel('Power (nA^2/Hz)')
        xlabel('Frequency (Hz)')
        hold all
    end

    % bring figure to front
    figure(hf);
    hax = get(hf,'Children');
    freqs = freqs(2:imax)';
    dfft = dfft(2:imax,:);
    % and plot
    plot(hax,freqs,dfft);
end

