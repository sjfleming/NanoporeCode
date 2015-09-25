% Stephen Fleming
% 4/13/15
% Quickly analyzes low-voltage diffusion experiment data, thoroughly.

function [blocks, open_pore] = blockage_analysis(sigdata,tr,cr)
    %% analyze
    
    % parameters
    %captureVoltage = 160; % in mV
    %holdingVoltageRange = [20 90]; % in mV
    %conductanceCutoff = 1.5; % in nS, if captured molecule conductance is greater than this, it's thrown out
    
    warning('off','signal:findpeaks:largeMinPeakHeight') % supress warning output
    display('Finding likely event starts...')
    
    % low-res data
    data = sigdata.getViewData(tr); % get downsampled data
    voltageCutoff = 10; % minimum voltage transition dV
    dV = [diff(data(:,3)); 0];
    
    % find all possible capture events
    [~,end_view_ind] = findpeaks(-1*dV.*(dV<(-1*voltageCutoff)),'MinPeakDist',5);
    potentialEnds = data(end_view_ind,1); % end times of potential events
    
    % eliminate ones that aren't real
    postVoltage = arrayfun(@(x) mean(data(x+20:x+30,3)), end_view_ind);
    duration = min([end_view_ind(1)-1, end_view_ind(1) - find(data(1:end_view_ind(1),3)>172,1,'last')]);
    duringVoltagePositive = arrayfun(@(x) prod(data((x-duration):x,3)>10), end_view_ind);
    logic = ( round(postVoltage/10)*10 == -10 & duringVoltagePositive );
    potentialEnds = potentialEnds(logic);
    
    %% go through and find exact information
    blocks = cell(1,numel(potentialEnds));
    open_pore = cell(1,numel(potentialEnds));
    timer = nan(1,numel(potentialEnds));
    display('Collating information from each event...')
    
    for i = 1:numel(potentialEnds)
        tic;
        
        % go back to the beginning (V>172)
        end_ind = potentialEnds(i)/sigdata.si;
        start_ind = sigdata.findPrev(@(x) x(:,3)>171, end_ind);
        
        % find each level
        dVregion = diff(sigdata.get(round(start_ind-(0.005/sigdata.si)):round(end_ind+(0.005/sigdata.si)),3));
        [~,step_ind] = findpeaks(-1*dVregion.*(dVregion<(-0.5)),'MinPeakDist',(0.005/sigdata.si)); % indices of voltage step transitions
        step_ind = step_ind + round(start_ind-(0.005/sigdata.si));
        step_duration = median(diff(step_ind)); % in number of points
        
        if step_duration*sigdata.si<0.15 % otherwise something went wrong
            
            % measure current for middle 50% of each level, or the piece of
            % the level before molecule jumps out
            mean_current = nan(numel(step_ind),1);
            std_current = nan(numel(step_ind),1);
            voltage = nan(numel(step_ind),1);
            open_mean_current = nan(numel(step_ind),1);
            open_std_current = nan(numel(step_ind),1);
            open_voltage = nan(numel(step_ind),1);
            for j = 1:numel(step_ind)
                conductance = sigdata.get((step_ind(j)-step_duration):step_ind(j),2)./ ...
                    sigdata.get((step_ind(j)-step_duration):step_ind(j),3) * 1000; % conductance whole level in nS
%                 figure(2)
%                 clf(2)
%                 plot(conductance)
%                 ylim([0 5])
                if and( mean(conductance(3:13))>cr(1), mean(conductance(3:13))<cr(2) ) % if start of level is one blockage
                    if and( mean(conductance((end-13):(end-3)))>cr(1), mean(conductance((end-13):(end-3)))<cr(2) ) % if end of level is one blockage
                        % mean of middle 80% of level
                        ind1 = step_ind(j)-round(0.9*step_duration);
                        ind2 = step_ind(j)-round(0.1*step_duration);
                    else
                        % find the real end of level and do middle 80% of
                        % what's left
                        this_step_ind = sigdata.findPrev(@(x) (x(:,2)*1000./x(:,3))<cr(2) & (x(:,2)*1000./x(:,3))>cr(1), step_ind(j));
                        this_step_duration = step_duration - (step_ind(j)-this_step_ind);
                        ind1 = this_step_ind-round(0.9*this_step_duration);
                        ind2 = this_step_ind-round(0.1*this_step_duration);
                    end
                    mean_current(j) = mean(sigdata.get(ind1:ind2,2));
                    std_current(j) = std(sigdata.get(ind1:ind2,2));
                    voltage(j) = mean(sigdata.get(ind1:ind2,3));
%                     hold on
%                     line([ind1 ind1]-(step_ind(j)-step_duration),[0 5],'Color','k')
%                     line([ind2 ind2]-(step_ind(j)-step_duration),[0 5],'Color','k')
%                     plot(mean([ind1 ind2]-(step_ind(j)-step_duration)),mean_current(j)*1000/voltage(j),'or','MarkerSize',10)
                elseif and( mean(conductance(3:13))>cr(3), mean(conductance((end-13):(end-3)))>cr(3) ) % check if we have an open pore level
                    % save this data for open pore
                    ind1 = step_ind(j)-round(0.9*step_duration);
                    ind2 = step_ind(j)-round(0.1*step_duration);
                    
                    open_mean_current(j) = mean(sigdata.get(ind1:ind2,2));
                    open_std_current(j) = std(sigdata.get(ind1:ind2,2));
                    open_voltage(j) = mean(sigdata.get(ind1:ind2,3));
                end
%                 pause()
            end
            
        end
        
        % save stuff
        blocks{i}.start = sigdata.get(start_ind,1);
        blocks{i}.finish = sigdata.get(end_ind,1);
        blocks{i}.start_ind = start_ind;
        blocks{i}.finish_ind = end_ind;
        blocks{i}.mean_current = mean_current;
        blocks{i}.std_current = std_current;
        blocks{i}.voltage = voltage;
        
        open_pore{i}.start = sigdata.get(start_ind,1);
        open_pore{i}.finish = sigdata.get(end_ind,1);
        open_pore{i}.start_ind = start_ind;
        open_pore{i}.finish_ind = end_ind;
        open_pore{i}.mean_current = open_mean_current;
        open_pore{i}.std_current = open_std_current;
        open_pore{i}.voltage = open_voltage;
        
        timer(i) = toc;
        % display progress
        if mod(i,10)==0
            display([num2str(round((numel(potentialEnds)-i)*nanmean(timer))) ' seconds left'])
        end
    end
    display('Done.')
    
end

