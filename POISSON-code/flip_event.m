function event = flip_event( event )

    event.mean = flipud(event.mean);
    event.stdv = flipud(event.stdv);
    event.length = flipud(event.length);
    event.start = flipud(event.start);
    event.ref_align = flipud(event.ref_align);
    
    % flip the model too
    event.model = flip_model(event.model);
    
    % and the sequence for realz
    event.sequence = seqrcomplement(event.sequence);

end

