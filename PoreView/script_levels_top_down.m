%% Stephen Fleming
% 6/23/15

% Top-down level extraction procedure

%% load data
tic;
% tr = pv.getCursors();
tr = [900 1000];
data = pv.data.getByTime(tr);
raw = pv.data.getViewData(tr);
time = data(:,1);
current = data(:,2)*1000;
raw(:,2) = raw(:,2)*1000;
d = [time, current];
clear data

%% filter heavily and downsample

medianRepeats = 5;
targetFs = 200; % Hz, final frequency after downsampling
lp = filt_lpb(d,4,1000);
[filtered] = filt_severe_median_downsample(lp,medianRepeats,(1/pv.data.si)/targetFs);
time = filtered(:,1);
current = filtered(:,2);
clear d filtered

%% find level changes

% take derivative
smoothing = 15;
deriv = smooth_derivative(current, smoothing);

% plot
figure(4)
clf(4)
plot(raw(:,1), raw(:,2))
hold on
plot(time, current,'LineWidth',2)
plot(time, deriv)
ylabel('Current (pA)')
xlabel('Time (s)')
legend('Raw data','Filtered data','Derivative')
xlim([min(time) max(time)])
set(gca,'FontSize',18)

%%

level_currents = mean(current);
level_timing = [time(1) time(end)];
pValue = 1e-5;

% find peaks in the derivative and sort them
[deriv_peaks, deriv_peak_locs] = findpeaks(abs(deriv),'MinPeakDist',2*smoothing,'MinPeakHeight',0.25*std(deriv));
[~, sort_ind] = sort(deriv_peaks,'descend');
num = 1;

%while metric(end)>0.15
for i = 1:numel(deriv_peaks)
    
    % find where level would split
    split_time = time(deriv_peak_locs(sort_ind(num)));
    split_level_ind = find(split_time>level_timing(:,1),1,'last');
    
    % check to make sure levels are statistically different
    ts_ind1 = find(time>level_timing(split_level_ind,1),1,'first');
    tf_ind1 = find(time<level_timing(split_level_ind,2),1,'last');
    [~,p,stat] = kstest2(current(ts_ind1:deriv_peak_locs(sort_ind(num))),current(deriv_peak_locs(sort_ind(num)):tf_ind1));
    %display(num2str(stat))
    
    if p<pValue
        
        % split off levels
        if split_level_ind==1
            level_timing = [level_timing(1,:); level_timing];
            level_currents = [level_currents(1); level_currents];
        elseif split_level_ind==numel(level_currents)
            level_timing = [level_timing; level_timing(end,:)];
            level_currents = [level_currents; level_currents(end)];
        else
            level_timing = [level_timing(1:split_level_ind,:); level_timing(split_level_ind,:); level_timing((split_level_ind+1):end,:)];
            level_currents = [level_currents(1:split_level_ind); level_currents(split_level_ind); level_currents((split_level_ind+1):end)];
        end
        level_timing(split_level_ind,2) = split_time;
        level_timing((split_level_ind+1),1) = split_time;
        
        % update current of new level 1
        ts_ind = find(time>level_timing(split_level_ind,1),1,'first');
        tf_ind = find(time<level_timing(split_level_ind,2),1,'last');
        level_currents(split_level_ind) = mean(current(ts_ind:tf_ind));
        
        % update current of new level 2
        ts_ind2 = find(time>level_timing(split_level_ind+1,1),1,'first');
        tf_ind2 = find(time<level_timing(split_level_ind+1,2),1,'last');
        level_currents(split_level_ind+1) = mean(current(ts_ind2:tf_ind2));
        
        discrete_current = ones(numel(time),1)*level_currents(1);
        for j = 2:numel(level_currents)
            ind = find(time>level_timing(j,1));
            discrete_current(ind:end,1) = level_currents(j);
        end
        
    end
    
    % plot levels
    figure(5)
    clf(5)
    plot(raw(:,1),raw(:,2),'k')
    hold on
    %line(level_timing',(ones(2,1)*level_currents'),'LineWidth',2)
    plot(time,discrete_current,'LineWidth',2)
    ylabel('Current (pA)')
    xlabel('Time (s)')
    xlim([min(time) max(time)])
    ylim([30 90])
    set(gca,'FontSize',18)
    %pause();
    metric(i) = sum((discrete_current - current).^2)/numel(time);
    %display(num2str(i))
    num = num+1;
    
end

%%
figure(6)
clf(6)
numlevels = 2:(numel(metric)+1);
plot(numlevels,metric,'o-')
hold on
plot(numlevels,metric+sqrt(numlevels),'o-')
plot(numlevels,metric+sqrt(numlevels)/5,'o-')
title('Metrics for Determining Final Number of Levels')
ylabel('Metric')
xlabel('Number of Levels')
set(gca,'FontSize',18)

%[~,optimumLevels] = min(metric+sqrt(numlevels));
[~,optimumLevels] = find(metric<0.1,1,'first');
optimumLevels = optimumLevels+1;

display(['Optimum number of levels is ' num2str(optimumLevels)])
%%
% reset and redo to the optimal number of levels
level_currents = mean(current);
level_timing = [time(1) time(end)];
num = 1;
for i = 1:(optimumLevels-1)
    
    % find where level would split
    split_time = time(deriv_peak_locs(sort_ind(num)));
    split_level_ind = find(split_time>level_timing(:,1),1,'last');
    
    % check to make sure levels are statistically different
    ts_ind1 = find(time>level_timing(split_level_ind,1),1,'first');
    tf_ind1 = find(time<level_timing(split_level_ind,2),1,'last');
    [~,p,stat] = kstest2(current(ts_ind1:deriv_peak_locs(sort_ind(num))),current(deriv_peak_locs(sort_ind(num)):tf_ind1));
    %display(num2str(stat))
    
    if p<pValue
        
        % split off levels
        if split_level_ind==1
            level_timing = [level_timing(1,:); level_timing];
            level_currents = [level_currents(1); level_currents];
        elseif split_level_ind==numel(level_currents)
            level_timing = [level_timing; level_timing(end,:)];
            level_currents = [level_currents; level_currents(end)];
        else
            level_timing = [level_timing(1:split_level_ind,:); level_timing(split_level_ind,:); level_timing((split_level_ind+1):end,:)];
            level_currents = [level_currents(1:split_level_ind); level_currents(split_level_ind); level_currents((split_level_ind+1):end)];
        end
        level_timing(split_level_ind,2) = split_time;
        level_timing((split_level_ind+1),1) = split_time;
        
        % update current of new level 1
        ts_ind = find(time>level_timing(split_level_ind,1),1,'first');
        tf_ind = find(time<level_timing(split_level_ind,2),1,'last');
        level_currents(split_level_ind) = mean(current(ts_ind:tf_ind));
        
        % update current of new level 2
        ts_ind2 = find(time>level_timing(split_level_ind+1,1),1,'first');
        tf_ind2 = find(time<level_timing(split_level_ind+1,2),1,'last');
        level_currents(split_level_ind+1) = mean(current(ts_ind2:tf_ind2));
%         display([ 'Level from [' num2str(time(ts_ind)) ' ' num2str(time(tf_ind)) '] mean is ' num2str(level_currents(split_level_ind)) ])
%         pause();
        
        discrete_current = ones(numel(time),1)*level_currents(1);
        for j = 2:numel(level_currents)
            ind = find(time>level_timing(j,1));
            discrete_current(ind:end,1) = level_currents(j);
        end
        
    end
    
    % plot levels
    figure(5)
    clf(5)
    plot(raw(:,1),raw(:,2),'k')
    hold on
    line(level_timing',(ones(2,1)*level_currents'),'LineWidth',2)
    %plot(time,discrete_current,'LineWidth',2)
    ylabel('Current (pA)')
    xlabel('Time (s)')
    xlim([min(time) max(time)])
    ylim([30 100])
    set(gca,'FontSize',18)
    %pause();
    metric(i) = sum((discrete_current - current).^2)/numel(time);
    %display(num2str(i))
    num = num+1;
    
end

title([num2str(numel(level_currents)) ' levels'])

figure(2)
clf(2)
plot(level_currents,'o-')
ylabel('Level Mean Current (pA)')
xlabel('Level Number')
set(gca,'FontSize',18)
title('Squiggle data')

toc;
