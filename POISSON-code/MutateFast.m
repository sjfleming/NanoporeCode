function [seq, events, nbases] = MutateFast(seq, seqs, events, params)

    nmut = numel(seqs);
    
    seqals = cell([nmut,1]);
    reflikes = cell([nmut,1]);
    dlikes = cell([nmut,1]);
    alscores = zeros([nmut,1]);
    
    % compare all mutations to best path, using fuller alignment
    params1 = params;
    params1.stripe_width = 400;
    tic
    [~,events,seqreflike] = align_likes(seq,events,params1);
    fprintf('Initial align took %0.1f\n',toc);

    fprintf('Aligning with best path');
    tic

    for i=1:numel(seqs)
        [alscore,seqals{i}] = seqalign(seq,seqs{i});
        %seqals{i} = fillinds(swfast(seq,seqs{i},[1 numel(seq);1 numel(seqs{i})],1000));
        % display alignment score and how many were aligned
        fprintf('.');
        fprintf('%0.1f|%d,',alscore,seqals{i}(end,2)-seqals{i}(1,2));
        alev = mapaligns(events,seqals{i});
        [scores,~,reflikes{i}] = align_likes(seqs{i},alev,params);
        alscores(i) = sum(scores);
        dlike = seqals{i};
        dlike = dlike - 1;
        dlike = dlike(and(dlike(:,1)>0,dlike(:,2)>0),:);
        %{
        dreflik = [0; diff(seqreflike)];
        dlike(:,1) = dreflik(dlike(:,1));
        dreflik = [0; diff(reflikes{i})];
        dlike(:,2) = dreflik(dlike(:,2));
        %}
        dlike(:,1) = seqreflike(dlike(:,1));
        dlike(:,2) = reflikes{i}(dlike(:,2));
        dlike = [0 0; diff(dlike)];
        dlikes{i} = cusum(dlike(:,2) - dlike(:,1));
        % remove areas where the likelihoods are the same
        % since these are probably identically-aligned regions
        sames = abs(dlike(:,2)-dlike(:,1))<1e-5;
        sames = imerode(sames,[1;1;1]) > 0.5;
        dlikes{i}(sames) = 0;
    end
    %fprintf('\b)\n');
    fprintf('\n');
    fprintf('Alignments took %0.1f\n',toc);
       
    % now let's just start mutating, or something
    mutations = [];
    
    tic
    ntot = 0;
    while numel(mutations) < numel(seq)/3
        dlms = zeros(numel(dlikes),1);
        % pick strand with largest dlikes
        for j=1:numel(dlikes)
            dlms(j) = max(dlikes{j});
        end
        [~,i] = max(dlms);

        dlike = dlikes{i};
        % get max index of dlikes
        [m,ind] = max(dlike);
        if m < 0.5
            % we're definitely done here
            break
        end
        % and prev and next zero
        i0 = find(and(dlike'==0,1:numel(dlike)<ind),1,'last');
        i1 = find(and(dlike'==0,1:numel(dlike)>ind),1,'first');
        if isempty(i0)
            i0 = 1;
        end
        if isempty(i1)
            i1 = numel(dlike);
        end
        % convert them to indices into both strands
        inds = seqals{i}([i0 ind i1],:);
        sind0 = inds(1,1);
        sind1 = inds(2,1);
        % get original and mutated seq fragment, and tighten them up
        mutpiece = seqs{i}(inds(1,2):inds(2,2));
        origpiece = seq(sind0:sind1);
        % remove identical bases and trim them down
        while numel(mutpiece) > 0 && numel(origpiece) > 0 && mutpiece(1) == origpiece(1)
            mutpiece = mutpiece(2:end);
            origpiece = origpiece(2:end);
            sind0 = sind0 + 1;
        end
        % from the other direction too
        while numel(mutpiece) > 0 && numel(origpiece) > 0 && mutpiece(end) == origpiece(end)
            mutpiece = mutpiece(1:end-1);
            origpiece = origpiece(1:end-1);
            sind1 = sind1 - 1;
        end
        % if we didn't find any, or if already mutated, or they're the
        % same, don't do nothin
        if ~isempty(sind0) && ~isempty(sind1) && ~strcmp(mutpiece,origpiece) %&& (numel(origpiece) > 1 || numel(mutpiece) > 1)
            % create our little mutation struct
            mut = [];
            mut.start = sind0;
            mut.original = origpiece;
            mut.mutation = mutpiece;
            if isempty(mutations)
                mutations = mut;
            else
                ntot = ntot + 1;
                mutations(end+1) = mut;
            end
        end
        % lastly, set that range of dlikes to 0 so we don't use it over
        dlikes{i}(i0:i1) = 0;
    end
    
    tic
    [seq,events,nbases] = MutateSequence(seq, events, params, mutations);
    fprintf('Mutations took %0.1f with %d\n',toc,nbases);

    % realign, but not full-frame
    
    fprintf('\n');
end