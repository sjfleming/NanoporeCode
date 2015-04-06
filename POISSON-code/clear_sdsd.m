function events = clear_sdsd(events)
% Clears standard deviation info from event models

    for i=1:numel(events)
        foo = 0*events(i).model.sd_stdv;
        events(i).model.sd_stdv = foo + 1/sqrt(2*pi);
        events(i).model.sd_mean = foo + 0;
        events(i).stdv = 0*events(i).stdv;
    end

end