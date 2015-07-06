dirbase = 'C:\Minion\Laszlo\Data\';
datadirs = dir([dirbase '\*phixfull*']);
datadirs = {datadirs.name};
for i=1:numel(datadirs)
    datadirs{i} = [dirbase datadirs{i}];
end
%%
% go through and plot or something
for n=1:numel(datadirs)

    load([datadirs{n} '\jumps_ian.mat'])
    load([datadirs{n} '\reduced.mat'])

    % plot the reduced data
    clf
    plot(reduced.pt/reduced.fSamp,reduced.data);
    hold on
    % and the jumps or something
    for i=1:numel(jumps.start)
        levels = jumps.ianjump{i}.means;
        %durs = 10*double(jump.durs+2)/reduced.fSamp;
        ys = doublemat(levels');
        %ts = doublemat([0 cumsum(durs)]');
        ts = jumps.reducedStart{i};
        ts(end+1) = ts(end) + jumps.duration{i}(end);
        ts = doublemat(rshape(ts));
        ts = ts(2:end-1);
        ts = reduced.pt(ts) / reduced.fSamp;
        plot(ts,ys,'r')
    end
    pause
end
%%
% or go through and save events
events = [];
for n=1:numel(datadirs)

    try
        load([datadirs{n} '\jumps_ian.mat'])
    catch
        continue
    end

    % and the jumps or something
    for i=1:numel(jumps.start)
        % make the event n stuff
        event = [];
        event.mean = jumps.ianjump{i}.means;
        event.stdv = jumps.std{i};
        event.median = jumps.median{i};
        event.start = jumps.reducedStart{i};
        event.count = numel(event.mean);
        if isempty(events)
            events = event;
        else
            events(end+1) = event;
        end
    end
    fprintf('%d: %d\n',n,numel(events));
end