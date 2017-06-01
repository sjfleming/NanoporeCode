function [xx,yy] = plot_hist(sigdata, tr, chan)
% plot a histogram of data in the channel specified by chan, from sigdata,
% in the time range tr

    % is this current or conductance?
    signames = sigdata.getSignalList();
    isCond = false;
    if ~isempty(strfind(signames{chan},'Conductance')) % if the word 'Conductance' is in the signal name
        isCond = true;
    end
    
    % create the x-values
    rough_downsampled_data = sigdata.getViewData(tr);
    rough_downsampled_sig = medfilt1(rough_downsampled_data(:,chan)*(1000-999*double(isCond)),50); % x1000 to get pA, if it's not conductance
    alpha = 10; % for recording from axopatch
    digitization = 1/alpha * 0.30517578125; % in pA
    if isCond
        maxVoltage = mode(abs(round(medfilt1(rough_downsampled_data(:,3),50)))); % most common voltage
        digitization = digitization/maxVoltage;
    end
    xx = min(rough_downsampled_sig)-500*digitization:digitization:max(rough_downsampled_sig)+500*digitization;
    
    % figure out how many chunks of data we should divide up the set into
    chunk_size = 2^22; % reasonable # of points
    num_chunks = ceil(diff(tr)/sigdata.si/chunk_size);

    % loop through and grab data
    yy = zeros(size(xx));
    ind = tr(1)/sigdata.si + 1; % absolute file index
    maxpt = tr(2)/sigdata.si + 1; % the last point to grab
    figure()
    for i = 1:num_chunks

        % grab that chunk's data
        d = sigdata.get(ind:min(maxpt,ind+chunk_size));
        data = d(d(:,chan)~=0,chan)*(1000-999*double(isCond)); % x1000 to get pA, if it's not conductance

        % add that binned histogram data to a running total
        yy_for_this_chunk = hist(data,xx);
        yy = yy + yy_for_this_chunk; % add to cumulative histogram

        % plot intermediate steps
        bar(xx,yy)
        drawnow;

        ind = ind + chunk_size + 1; % increment for next time
    end

    % finish the figure
    bar(xx,yy) % bar graph
    ylabel('Number of data points')
    if isCond
        xlabel('Conductance (nS)')
    else
        xlabel('Current (pA)')
    end
    title([sigdata.filename(end-18:end-4) ' [' num2str(round(tr(1)*10)/10) ...
        ',' num2str(round(tr(2)*10)/10) '], ' signames{chan}],'interpreter','none')
    set(gca,'fontsize',18,'looseinset',[0 0 0 0],'outerposition',[0.01 0.01 0.98 0.98])

end