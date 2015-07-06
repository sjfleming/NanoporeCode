% Stephen Fleming
% 4/13/15
% Quickly analyzes low-voltage diffusion experiment data, thoroughly.

function blocks = blockage_analysis(sigdata,tr)
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
    potentialEnds = data(end_view_ind,1); % start times of potential events
    
    % eliminate ones that aren't real
    postVoltage = arrayfun(@(x) mean(data(x+20:x+30,3)), end_view_ind);
    logic = ( round(postVoltage/10)*10 == -10 );
    potentialEnds = potentialEnds(logic);
    
    %% go through and find exact information
    blocks = cell(1,numel(potentialEnds));
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
            
            % measure current for middle 50% of each level
            mean_current = nan(numel(step_ind),1);
            voltage = nan(numel(step_ind),1);
            for j = 1:numel(step_ind)
                mean_current(j) = mean(sigdata.get((step_ind(j)-round(0.75*step_duration)):(step_ind(j)-round(0.25*step_duration)),2));
                voltage(j) = mean(sigdata.get((step_ind(j)-round(0.75*step_duration)):(step_ind(j)-round(0.25*step_duration)),3));
            end
            
        end
        
        % save stuff
        blocks{i}.start = sigdata.get(start_ind,1);
        blocks{i}.finish = sigdata.get(end_ind,1);
        blocks{i}.start_ind = start_ind;
        blocks{i}.finish_ind = end_ind;
        blocks{i}.mean_current = mean_current;
        blocks{i}.voltage = voltage;
        
        timer(i) = toc;
        % display progress
        if mod(i,10)==0
            display([num2str(round((numel(potentialEnds)-i)*nanmean(timer))) ' seconds left'])
        end
    end
    display('Done.')
    
end

