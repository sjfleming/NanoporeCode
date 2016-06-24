function [seq, events] = MutateExtend(seq, events, params)

    % compare all mutations to best path, using fuller alignment
    params1 = params;
    params1.stripe_width = 400;
    tic
    [~,events] = align_likes(seq,events,params1);
    fprintf('Initial align took %0.1f\n',toc);
    
    mutations = [];

    nbases = 4;
    
    % create all mutations of length 4 (256 of them)
    for i=1:4^nbases
        lbits = bitget(i-1,1:2:2*nbases);
        hbits = bitget(i-1,2:2:2*nbases);
        bases = int2nt(1+2*hbits+lbits);
        
        mut = [];
        mut.start = numel(seq);
        mut.original = [];
        mut.mutation = bases;
        if isempty(mutations)
            mutations = mut;
        else
            mutations(end+1) = mut;
        end
        mut = [];
        mut.start = 1;
        mut.original = [];
        mut.mutation = bases;
        mutations(end+1) = mut;
    end

    tic
    [seq, events, nbases] = MutateSequence(seq, events, params, mutations);
    fprintf('Mutation testing took %0.1f with %d\n',toc,nbases);
end