function filtevents = filter_events(seq, events)

    filtevents = [];
    for i=1:numel(events)
        [score, dpath] = seqalign(seq,events(i).sequence);
        if score < 60 || size(dpath,1) < 500
            continue
        end
        if isempty(filtevents)
            filtevents = events(i);
        else
            filtevents(end+1) = events(i);
        end
    end
    
    fprintf('Filtered out %d events\n',numel(events)-numel(filtevents));

end