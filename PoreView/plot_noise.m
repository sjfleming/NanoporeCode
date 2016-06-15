function plot_noise(sigdata, trange)
%PLOT_NOISE Makes or adds to existing noise plot
%   plot_noise(sigdata, trange)
%   Takes a SignalData object and a time range as input, uses the same
%       algorithm as ClampFit (I think).

    % do same thing as ClampFit does - average spectral segs
    % ClampFit does 2*65536 pts per seg (2^17)
    % get start and end index
    irange = floor(trange/sigdata.si);
    
    fftsize = min(diff(irange),2*65536);

    % only process real signals
    dfft = zeros(fftsize,sigdata.nsigs);
    % number of frames
    nframes = 0;
    
    wh = waitbar(0,'Calculating power spectrum...','Name','PoreView');
    
    for ind=irange(1):fftsize:irange(2)
        % get only the real signals
        d = sigdata.get(ind:min(irange(2),ind+fftsize-1),1+(1:sigdata.nsigs));
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
    %temp = input('Temperature in degrees C: '); % in degrees C
    temp = 25;

    % do the averaging
    dfft = dfft / nframes;
    
    % dave's factor of 2
    dfft = 2*dfft;
    
    % and calculate the frequency range
    f = 1/sigdata.si*(0:fftsize-1)/fftsize;
    
    % range of f to plot, only do half (Nyquist and all)
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
        set(ax,'XLim',[1 f(imax)]);
        set(gca,'FontSize',20)
        set(gca,'LooseInset',[0 0 0 0]) % the all-important elimination of whitespace!
        set(gca,'OuterPosition',[0 0 0.99 1]) % fit everything in there
        set(gcf,'Position',[100 500 750 500]) % size the figure
        % grid
        grid on
        grid minor
        box on
        % labels
        title('Noise Power Spectrum');
        ylabel('Current Noise (nA^2/Hz)')
        xlabel('Frequency (Hz)')
    end

    % bring figure to front
    figure(hf);
    % and plot
    V = mean(sigdata.get(trange(1)/sigdata.si:(trange(1)/sigdata.si+1000),3));
    I = mean(sigdata.get(trange(1)/sigdata.si:(trange(1)/sigdata.si+1000),2))*1000;
    if abs(V)<4
        conductance = input('Conductance (nS): ');
    else
        conductance = I/V; % in nS
        %display(['Conductance = ' num2str(conductance,3) 'nS'])
    end

    R = 1/(conductance*1e-9);
    k = 1.38 * 10^-23; % Boltzmann constant
    T = temp + 273.15; % absolute temperature in Kelvin
    Rf = 500e6; % feedback resistor in Axopatch with beta = 1 in whole cell mode
    Cin = 5.5e-12; % headstage input capacitance is about 4pF
    V_headstage = 10e-9; % input voltage noise on headstage op-amp = 3nV/sqrt(Hz)
    Ra = 5.6e6; % access resistance = rho/4*a (J.E. Hall, 1975), rho = 0.0895ohm*m for 1M KCl (http://www.sigmaaldrich.com/catalog/product/fluka/60131?lang=en&region=US), a = 1nm
    Cm = 0.25e-12; % typical membrane capacitance 0.2pF
    johnson = 4*k*T * ( 1/R + 1/Rf ) * 10^18;
    shot = 0; %2 * 1.602*1e-19 * abs(I) * 1e-12 * 10^18;
    ion = 1.12e-9 * (I*1e-12)^2 * (1e9)^2; % ion number fluctuations proportional to current squared (10/20/16 data for MspA from Oxford)
    Y = @(f,Ra,Cm) sqrt(-1)*2*pi*f*Cin + 1./( Ra + 1./( 1/R + sqrt(-1)*2*pi*f.*Cm ) ) + 1./Rf; % complex conductance (admittance)
    headstage = V_headstage^2*real(Y(f,Ra,Cm).*conj(Y(f,Ra,Cm))) * 10^18;
    % the filter
    cutoff = 1e4;
    factor = 2/1.8031;
    [b,a] = besself(4,2*cutoff/factor);
    h = freqs(b,a,f);
    noise_model = 4*k*T*real(Y(f,Ra,Cm)*1e18.*h.*conj(h)) + ion*real(h.*conj(h)) + shot*real(h.*conj(h)) + headstage.*real(h.*conj(h));
    %display(' ')
    %display(['V = ' num2str(V,4) ' mV'])
    %display(['I = ' num2str(I,4) ' pA'])
    %display(['Johnson noise = ' num2str(johnson,3) ' nA^2/Hz'])
    %display(['Shot noise = ' num2str(shot,3) ' nA^2/Hz'])
    %display(['Sum = ' num2str(johnson+shot,3) ' nA^2/Hz'])
    i1 = find(f>80,1,'first');
    i2 = find(f>280,1,'first');
    measured_white = mean(dfft(i1:i2,1));
    measured_white_err = std(dfft(i1:i2,1))/sqrt(i2-i1+1);
    %display(['Measured white noise = ' num2str(measured_white,3) ' ± ' num2str(measured_white_err,2) ' nA^2/Hz'])
    display([num2str(V,4) ', ' num2str(I,4) ', ' num2str(measured_white,3) ', ' num2str(measured_white_err,2) ';'])
    rms = sqrt(mean(dfft(i1:i2,1))*1e6*(f(i2)-f(i1)));
    rms_dev = sqrt(std(dfft(i1:i2,1))*1e6*(f(i2)-f(i1))/sqrt(i2-i1+1));
    %display(['Irms = ' num2str(rms,3) ' ± ' num2str(rms_dev,2) ' pA over a ' num2str(round(f(i2)-f(i1))) 'Hz bandwidth'])
    line([1 1e5],johnson*ones(1,2),'Color','k','LineStyle','--')
    hold on
    plot(f(1:imax)',dfft(1:imax,1));
    line(f(1:imax),noise_model(1:imax),'Color','c','LineStyle','-')
    name = [sigdata.filename(65:68) '\_' sigdata.filename(70:71) '\_' sigdata.filename(73:74) '\_' sigdata.filename(76:end-4)];
    annotation('textbox', [0.75 0.9 0 0], 'String', name, 'FontSize', 20);
    %ylim([0.2 2.5]*1e-10)
    %xlim([1 10000])
%     I_list = evalin('base','I');
%     I_list = [I_list, I];
%     S_list = evalin('base','S');
%     S_list = [S_list, measured_white];
%     err_list = evalin('base','err');
%     err_list = [err_list, measured_white_err];
%     assignin('base','I',I_list);
%     assignin('base','S',S_list);
%     assignin('base','err',err_list);
end

