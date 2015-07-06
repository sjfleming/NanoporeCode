function model = flip_model(model)
% flips a model around to its reverse complement order

    % AGGTC - 0010101101
    % GACCT - 1000010111
    
    nstates = 1024;
    
    inds = zeros(nstates,1);
    
    for i=1:nstates
        % flip order of bits in groups of 2, and invert them
        bits = ~bitget(i-1,[2 1 4 3 6 5 8 7 10 9]);
        % now put back into integer
        inds(i) = 1+bits*2.^(9:-1:0)';
    end

    model.level_mean = model.level_mean(inds);
    model.level_stdv = model.level_stdv(inds);
    model.sd_mean = model.sd_mean(inds);
    model.sd_stdv = model.sd_stdv(inds);
    model.kmer = model.kmer(:,inds);    
%    model.weight = model.weight(inds);
%    model.variant = model.variant(inds);
end