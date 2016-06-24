function [discreteData] = findLevels(data,fs)

% Stephen Fleming
% 6/19/14

% findLevels(data,fs) returns the levels and durations in the current
% modulation event passed in as 'data'. Implements 'levelFilter.m'.
% fs is the sampling rate. Level durations are returned as a matrix where
% the first column is the start time of a level and the second column is
% the end time of that level. 'discreteData' has the first column as time
% points and the second as current.

%% Filter heavily

display('Heavily filtering...')
[subsampled,f] = levelFilter(data,fs);
%timef = linspace(0,numel(data)*1/fs,numel(f)); % create a time vector

%% Itentify level transitions using a difference signal

display('Identifying level transitions...')
% Take the difference at each point between the average of the next n
% and the average of the previous n points
n = 10;
for i = 2:numel(f)-1
    before = mean(f(max(i-n,1):i-1));
    after = mean(f(i+1:min(i+n,numel(f))));
    
    diff(i) = after-before;
end

% % Identify peaks in this difference signal
% [c,index] = findpeaks(abs(diff).*(abs(diff)-2*std(diff)>0));
% 
% %% Figure out levels and timing
% 
%discreteData = [timef' f']; % time vector
% 
% a = 1;
% for i = 1:numel(index)
%     levels(i,1) = mean(f(a:index(i))); % level itself
%     discreteData(2,a:index(i)) = levels(i,1)*ones(index(i)-a+1,1); % vector data
%     timing(i,1) = timef(a); % start time
%     timing(i,2) = timef(index(i)); % end time
%     a = index(i)+1;
% end


%% Estimate error and iterate to convergence

% display('Optimizing levels...')
% % Initialize loop
% j = 2;
% s = 1; % initialize prefactor
% errorMetric = [10000 9000];
% err = [];
% % Iterate until change in error is less than a set tolerance
% %while errorMetric(j-1)-errorMetric(j)>0.01 && s>0.2
% while s>0.5
%     % Move peak detection limit
%     s = s-0.1;
%     % Re-find peaks in difference signal
%     [~,index] = findpeaks(abs(diff).*(abs(diff)-s*std(diff)>0));
%     % Re-calculate levels
%     a = 1;
%     for i = 1:numel(index)
%         levels(i,1) = mean(f(a:index(i))); % level itself
%         discreteData(2,a:index(i)) = levels(i,1)*ones(index(i)-a+1,1); % vector data
%         timing(i,1) = timef(a); % start time
%         timing(i,2) = timef(index(i)); % end time
%         err(i) = sum(discreteData(2,a:index(i)) - f(a:index(i)))*1e10; % avaraged error (no absolute value) in this level
%         a = index(i)+1;
%     end
%     % Calculate error
%     errorMetric(j+1) = sum(abs(err));
%     numLevels(j-1) = numel(index);
%     display(['Iteration ' num2str(j-1) '. Error = ' num2str(errorMetric(j+1))])
%     j = j+1;
% end
% 
% display('Done.')

[~,index] = findpeaks(abs(diff).*(abs(diff)-3*std(diff)>0),'minpeakdistance',1000);

%% Merging levels which are statistically the same

% Do Kolmogorov-Smirnov test
inds = [1 index numel(f)];

worsts = [numel(inds) 1];

%pValues = [inds(2:end-1)' zeros(numel(inds)-2,1)];
pValues = [inds' zeros(numel(inds),1)];
for i = 2:numel(inds)-1
    [~,pValues(i,2)] = kstest2(subsampled(inds(i-1):inds(i)),subsampled(inds(i):inds(i+1)));
end
% for i = 2:numel(inds)-1
%     [~,pValues(i-1,2)] = kstest2(subsampled(inds(i-1):inds(i)),subsampled(inds(i):inds(i+1)));
% end

% Iterate until the largest p value falls below cutoff
while worsts(end,2) > 0.001 && size(pValues,1) > 2
    [pBad,indsind] = max(pValues(:,2)); % find worst p-value
    worsts(end+1,:) = [numel(inds) pBad]; % keep track of it
    %inds = inds((1:numel(inds))~=indsind+1); % remove that index from the list of transitions
    inds(indsind) = [];
    pValues(indsind,:) = [];
%     pValues(indsind,:) = NaN;
%     pValues = nonzeros(~isnan(pValues).*pValues);  % remove it from p-values
    if indsind == 2
        first = indsind;
    else
        first = indsind-1;
    end
    if indsind == size(pValues,1)
        second = indsind-1;
    else
        second = indsind;
    end
    if first < second
        for i = first:second  % recalculate the ones before and after (if there is a before and after)
            %display(['i=' num2str(i) ', indsind=' num2str(indsind) ', size(pValues)=' num2str(size(pValues))])
            [~,pValues(i,2)] = kstest2(subsampled(inds(i-1):inds(i)),subsampled(inds(i):inds(i+1)));
        end
    end
    %display([num2str(worsts(end,2))])
end

% figure(2)
% plot(worsts(:,1),worsts(:,2))

%% Look at the remaining merged levels

for i = 1:numel(inds)-1
    discreteData(inds(i):inds(i+1)) = mean(subsampled(inds(i):inds(i+1)));
end

% figure(1)
% hold on
% plot(discreteData(:,1),discreteData(:,2),'y')

end
