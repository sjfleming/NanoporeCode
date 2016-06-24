function events = clearaligns(events)

    for i=1:numel(events)
        events(i).ref_align = 0*events(i).ref_align;
    end

end