function events = mapaligns(events, dpath)
% Maps ref_aligns from one sequence to another using a SW-aligned dpath

    % first, make first row of dpath uniques, so we can do interp with it
    [~,ia] = unique(dpath(:,1));
    dpath = dpath(ia,:);

    for i=1:numel(events)
        % now we need to find the inds of refal in dpath
        %{
        for j=1:numel(events(i).ref_align)
            z = find(dpath(:,1)==events(i).ref_align(j),1,'first');
            if isempty(z)
                events(i).ref_align(j) = -1;
            else
                events(i).ref_align(j) = dpath(z,2);
            end
        end
        %}
        refal = events(i).ref_align;
        ra0 = refal > 0;
        events(i).ref_align(:) = -1;
        % get the points using interp1 with nearest
        % extrapolating to 0 outside function
        events(i).ref_align(ra0) = round(interp1(dpath(:,1),dpath(:,2),refal(ra0),'nearest',0));
        % and set extras to zero
        i0 = find(events(i).ref_align > 0,1,'first');
        i1 = find(events(i).ref_align > 0,1,'last');
        events(i).ref_align(1:i0-1) = 0;
        events(i).ref_align(i1+1:end) = 0;
    end

end