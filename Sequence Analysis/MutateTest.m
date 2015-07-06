function seq = MutateTest(seq, events, params)

    % compare all mutations to best path, using fuller alignment
    params1 = params;
    params1.stripe_width = 400;
    tic
    [~,events] = align_likes(seq,events,params1);
    fprintf('Initial align took %0.1f\n',toc);
    
    mutations = [];

    % create our little mutation struct
    for i=1:500
        
        % in a random place
        istart = randi(numel(seq)-50);
        % with length of a few bases
        num = randi([1 8]);
        
        mut = [];
        mut.start = istart;
        mut.original = seq(istart:istart+num-1);
        mut.mutation = mut.original;
        if isempty(mutations)
            mutations = mut;
        else
            mutations(end+1) = mut;
        end
    end

    tic
    [seq,~,nbases] = MutateSequence(seq, events, params, mutations);
    fprintf('Mutation testing took %0.1f with %d\n',toc,nbases);
end