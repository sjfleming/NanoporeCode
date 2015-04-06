function events = order_events(seq, events)
% find the orders of the events

    fprintf('Ordering events');
    for j=1:numel(events)
        sfwd = seqalign(events(j).sequence(1:500),seq);
        srev = seqalign(seqrcomplement(events(j).sequence(1:500)),seq);
        fprintf('.');
        % is the strand probably backwards?
        if srev > sfwd
            events(j) = flip_event(events(j));
        end
    end
    fprintf('\n');

end