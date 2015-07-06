function out = remove_spikes(d,channel,trange,spike,indexV,isTemplate)
%NEXT_TRANSITION loops through data in one second chunks and finds the
%   time of the first transition to a new level
%   Stephen Fleming

    if isTemplate
        
        % generate spike template

        data = d.getByTime(trange); % grab a section with a near-ideal spike

        dvdt = data(2:end,3)-data(1:end-1,3); % voltage derivative
        ideal = abs((data(end,channel)-data(1,channel))/(data(end,3)-data(1,3)))*(data(:,3)-data(1,3))+data(1,channel); % idealized current after correction
        spike = data(:,channel)-ideal; % current spike due to voltage change

        figure(2)
        clf(2)
        plot(data(1:end-1,1),dvdt,'r')
        hold on
        plot(data(:,1),spike)
        title('dV/dt and current spike')
        xlabel('Time (s)')

        [~,indexV] = min(dvdt); % index of the voltage change timing
        
        % return
        out = [indexV; spike];
        
    else
        
        % correct signal using spike template
        
        newdvdt = d(2:end,2)-d(1:end-1,2); % voltage derivative
        %[delta,ind] = min(newdvdt); % index of the voltage change timing
        [~,ind] = findpeaks(-1*newdvdt,'MinPeakHeight',10);
        correction = zeros(numel(d(:,1)),1);
        
        for j=1:numel(ind) % add spikes at the right positions
            start = ind(j)-indexV+1;
%             display(['Correction  : ' num2str([max(1,start), min(start+numel(spike)-1,numel(d(:,1)))])])
%             display(['Spike  : ' num2str([max(1,-1*start+2), min(numel(spike),numel(d(:,1))-start+1)])])
            correction(max(1,start):min(start+numel(spike)-1,numel(d(:,1)))) = ...
                spike(max(1,-1*start+2):min(numel(spike),numel(d(:,1))-start+1)); % put the spikes in, not going past ends
        end
        
        % return
        out = d(:,channel)-correction;
        
    end

end
