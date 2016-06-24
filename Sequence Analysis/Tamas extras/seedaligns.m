function newevents = seedaligns(seq, events, params)
% Uses SW alignment to initialize the ref_aligns of all the events

    params.stripe_width = 400;
    
    fprintf('Seeding alignments');
    
    newevents = events(1:0);

    for i=1:numel(events)
        % if it has an alignment, do a wide realignment
        if max(events(i).ref_align) > 0
            [~,newevents(end+1)] = align_likes(seq,events(i),params);
            fprintf('.');
            continue;
        end
        % no alignment, no sequence... leave it out
        if isempty(events(i).sequence)
            fprintf('X');
            continue;
        end
        % else, seed the align by aligning to itself first
        events(i).ref_align = round(linspace(1,numel(events(i).sequence),numel(events(i).ref_align)))';
        [~,events(i)] = align_likes(events(i).sequence,events(i),params);
        % and then mapping the alignment
        [score,dpath] = seqalign(events(i).sequence,seq);
        if score > 60 && size(dpath,1) > 200
            % using realigned align_likes form
            [~,newevents(end+1)] = align_likes(seq,mapaligns(events(i),dpath),params);
            fprintf('.');
        else
            fprintf('X');
        end
    end
    fprintf('\n');

end