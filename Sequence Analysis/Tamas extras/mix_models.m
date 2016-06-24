function events = mix_models(events)

    models = events(1).model;
    if events(1).flipped
        models = flip_model(models);
    end
    modelnames = {events(1).model.name};
    modelcounts = 1;
    
    for i=2:numel(events)
        inds = strcmp(events(i).model.name,modelnames);
        if any(inds)
            modelcounts(inds) = modelcounts(inds) + 1;
        else
            if ~events(i).flipped
                models(end+1) = events(i).model;
            else
                models(end+1) = flip_model(events(i).model);
            end
            modelnames{end+1} = events(i).model.name;
            modelcounts(end+1) = 1;
        end
    end
    
    props = {'level_mean','level_stdv','sd_mean','sd_stdv'};
    
    % now make averaged models
    function newmodel=addmodels(mod1,mod2)
        for j=1:numel(props)
            p = props{j};
            newmodel.(p) = mod1.(p) + mod2.(p);
        end
    end

    mods = {[],[]};
    nmods = [0 0];

    for i=1:numel(modelnames)
        if ~isempty(strfind(modelnames{i},'template'))
            ind = 1;
        else
            ind = 2;
        end
        if isempty(mods{ind})
            mods{ind} = models(i);
        else
            mods{ind} = addmodels(mods{ind},models(i));
        end
        nmods(ind) = nmods(ind) + 1;
    end
    % and now reset the models
    for i=1:numel(events)
        if ~isempty(strfind(events(i).model.name,'template'))
            ind = 1;
        else
            ind = 2;
        end
        model = mods{ind};
        model.weight = events(i).model.weight;
        if events(i).flipped
            model = flip_model(model);
        end
        for j=1:numel(props)
            p = props{j};
            events(i).model.(p) = normlevels(model.(p),events(i).model.(p));
        end
    end
end