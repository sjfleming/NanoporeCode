% Stephen Fleming
% 4/13/15
% Quickly analyzes low-voltage diffusion experiment data, thoroughly.

function events = blockage_analysis(sigdata,tr,captureVoltage,holdingVoltageRange,conductanceCutoff,fcondSig)
    %% analyze
    
    % parameters
    %captureVoltage = 160; % in mV
    %holdingVoltageRange = [20 90]; % in mV
    %conductanceCutoff = 1.5; % in nS, if captured molecule conductance is greater than this, it's thrown out
    
    warning('off','signal:findpeaks:largeMinPeakHeight') % supress warning output
    display('Finding likely event starts...')
    
    % low-res data
    data = sigdata.getViewData(tr); % get downsampled data
    voltageCutoff = 20; % minimum voltage transition dV
    dV = [diff(data(:,3)); 0];
    
    % find all possible event starts
    [~,locs] = findpeaks(-1*dV.*(dV<(-1*voltageCutoff)),'MinPeakDist',5);
    potentialStarts = data(locs,1); % start times of potential events
    
    % eliminate ones with wrong capture voltage, wrong holding voltage,
    % wrong captured conductance
    preVoltage = arrayfun(@(x) mean(data(max(1,x-20):max(1,x-10),3)), locs);
    postVoltage = arrayfun(@(x) mean(data(x+20:x+30,3)), locs);
    capturedConductance = arrayfun(@(x) mean(data(max(1,x-5):max(1,x-2),2)*1000./data(max(1,x-5):max(1,x-2),3)), locs);
    logic = ( round(preVoltage/10)*10==captureVoltage & ...
        postVoltage>holdingVoltageRange(1) & postVoltage<holdingVoltageRange(2) & ...
        capturedConductance<conductanceCutoff );
    goodLocs = locs(logic);
    potentialStarts = potentialStarts(logic);
    
    %% go through and find exact information
    events = cell(1,numel(potentialStarts));
    display('Collating information from each event...')
    for i = 1:numel(potentialStarts)
        tic;
        
        % find exact start
        ind = sigdata.findNext(@(x) x(:,3)<(holdingVoltageRange(2)+10), round(potentialStarts(i)/sigdata.si)); % find where we're low
        ind = sigdata.findPrev(@(x) x(:,3)>(captureVoltage-10), ind); % then go back to when we were high
        startInd = ind;
        start = startInd * sigdata.si;
        
        % find exact end
        ind1 = sigdata.findNext(@(x) x(:,fcondSig(1))>conductanceCutoff, startInd + ceil(80e-6/sigdata.si)); % find end totally out
        devent = sigdata.get( (ind1-ceil(60e-6/sigdata.si)) : (ind1-ceil(30e-6/sigdata.si)), fcondSig(1) );
        cut = mean(devent) + 2*std(devent);
        ind = sigdata.findPrev(@(x) x(:,fcondSig(1))<cut, ind1); % backtrack to when current increase started
        finishInd = ind;
        finish = finishInd * sigdata.si;
        
        % load data
        cBefore = sigdata.get([round((start-0.009)/sigdata.si), startInd], fcondSig(1));
        cAfter = sigdata.get([round((start+1e-4)/sigdata.si), round(min(start+0.001,max(finish,start+1e-4))/sigdata.si)], fcondSig(1));
        Vhold = sigdata.get([(startInd+50),(finishInd+100)],3);
        Vcapture = sigdata.get([max(1,(startInd-100)),max(2,(startInd-50))],3);
        
        % save stuff
        events{i}.start = start;
        events{i}.finish = finish;
        events{i}.start_ind = startInd;
        events{i}.finish_ind = finishInd;
        events{i}.duration = (finish-start);
        events{i}.hold_conductance = mean(cAfter);
        events{i}.hold_cond_std = std(cAfter);
        events{i}.capture_conductance = mean(cBefore); % from 9ms before to 1ms before
        events{i}.capture_cond_std = std(cBefore);
        events{i}.hold_voltage = round(mean(Vhold(1:50))*10)/10;
        events{i}.capture_voltage = round(mean(Vcapture));
        events{i}.kicked_out = min(Vhold)<5;
        
        timer(i) = toc;
        % display progress
        if mod(i,10)==0
            display([num2str(round((numel(potentialStarts)-i)*mean(timer))) ' seconds left'])
        end
    end
    display('Done.')
    
end

