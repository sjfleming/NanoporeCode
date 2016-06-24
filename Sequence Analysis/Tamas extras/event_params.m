function events = event_params(events, params)
% set the skip and stay parameters of all the events to be their averages

    mskips = [0 0];
    mskips(1) = mean(arrayfun(@(x) x.model.skip_prob,events(1:2:end)));
    mskips(2) = mean(arrayfun(@(x) x.model.skip_prob,events(2:2:end)));
    mstays = [0 0];
    mstays(1) = mean(arrayfun(@(x) x.model.stay_prob,events(1:2:end)));
    mstays(2) = mean(arrayfun(@(x) x.model.stay_prob,events(2:2:end)));
    mextend = mstays;
    
    if nargin > 1
        mskips = [params.skip_t params.skip_c];
        mstays = [params.stay_t params.stay_c];
        mextend = [params.ext_t params.ext_c];
    end
    
    for j=1:2:numel(events)
        % template strands
        events(j).model.skip_prob = mskips(1);
        events(j).model.stay_prob = mstays(1);
        events(j).model.extend_prob = mextend(1);
        events(j+1).model.skip_prob = mskips(2);
        events(j+1).model.stay_prob = mstays(2);
        events(j+1).model.extend_prob = mextend(2);
    end

end