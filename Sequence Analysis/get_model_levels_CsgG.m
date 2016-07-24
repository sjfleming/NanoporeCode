function [pA, pA_std, std_mean] = get_model_levels_CsgG(seq, measured_levels)
% get_model_levels_CsgG(seq, measured_levels) returns current levels for
% sequence 'seq' that are predicted by Oxford Nanopore Technologies, Inc.
% for their CsgG pore
% scaled to mesh with the empirical CDF of the levels passed in
% Stephen Fleming, 7/23/16

    % Takes in a sequence in the form of 'ATCGTTCAAAGC...' with 'R'
    % denoting abasic sites

    % generate "states"
    states = get_states(nt2int(seq), 6); % k=6 for k-mer

    % get the current from Oxford's models
    d = load('CsgG_model.mat');
    pA = nan(size(states));
    std_mean = ones(size(states));
    % real parts
    states_real = states(imag(states)==0);
    pA(imag(states)==0) = d.CsgG_model.level_mean(states_real);
    pA_std(imag(states)==0) = d.CsgG_model.level_stdv(states_real);
    std_mean(imag(states)==0) = d.CsgG_model.sd_mean(states_real);
    
    % find levels not contained in the model
    inds = find( d.CsgG_model.weight(states_real) < 2 ); % states with weight less than 2 seem to be garbage levels
    % the (previous index divisible by 4) + 1 is the index of XXXXXA for the same 5-mer
    % average over the 5-mers of the same sequence that have non-negligible weights
    for j = inds % go through each model state
        i = states_real(j);
        rem = mod(i,4); % rem=1 means XXXXXA, rem=2 means XXXXXC, etc.
        weighted_avg_mean = sum(d.CsgG_model.level_mean((i-rem+1):(i-rem+4)) .* d.CsgG_model.weight((i-rem+1):(i-rem+4))) ...
            / sum(d.CsgG_model.weight((i-rem+1):(i-rem+4)));  % XXXXXA through XXXXXT
        weighted_avg_stdv = sum(d.CsgG_model.level_stdv((i-rem+1):(i-rem+4)) .* d.CsgG_model.weight((i-rem+1):(i-rem+4))) ...
            / sum(d.CsgG_model.weight((i-rem+1):(i-rem+4)));  % XXXXXA through XXXXXT
        real_indices = nonzeros((1:numel(states)).*(imag(states)==0));
        pA(real_indices(j)) = weighted_avg_mean;
        pA_std(real_indices(j)) = weighted_avg_stdv;
    end
    
    % for now, input the levels exactly from measured data (2016_01_21_0009 at 1185sec has all these levels)
    if sum(imag(states)~=0)>0 % if there are imaginary numbers, we have abasics at those spots
        % go through each spot and match it to a known current level
        for ind = nonzeros((1:numel(states)).*(imag(states)~=0))'
            kmer = seq(ind:ind+4);
            switch kmer
                case 'RRRRR'
                    pA(ind) = 105.7;
                    pA_std(ind) = 0.74;
                case 'RRRRT'
                    pA(ind) = 105.0;
                    pA_std(ind) = 0.64;
                case 'RRRTT'
                    pA(ind) = 94.4;
                    pA_std(ind) = 1.35;
                case 'RRTTT'
                    pA(ind) = 56.2;
                    pA_std(ind) = 5.79;
                case 'RTTTT'
                    pA(ind) = 37.8;
                    pA_std(ind) = 0.87;
                case 'GCGAR'
                    pA(ind) = 59.3;
                    pA_std(ind) = 1.06;
                case 'CGARR'
                    pA(ind) = 72.7;
                    pA_std(ind) = 1.56;
                case 'GARRR'
                    pA(ind) = 86.5;
                    pA_std(ind) = 1.12;
                case 'ARRRA'
                    pA(ind) = 93.0;
                    pA_std(ind) = 1.28;
                case 'RRRAC'
                    pA(ind) = 81.0;
                    pA_std(ind) = 2.42;
                case 'RRACA'
                    pA(ind) = 62.0;
                    pA_std(ind) = 2.97;
                case 'RACAT'
                    pA(ind) = 49.6;
                    pA_std(ind) = 1.20;
                otherwise
                    display('Unknown abasic k-mer in sequence: guessing!');
                    pA(ind) = (sum(kmer=='R')*105 + sum(kmer=='A')*53 + sum(kmer=='C')*52 + sum(kmer=='G')*55 + sum(kmer=='T')*36)/5;
                    pA_std(ind) = 1.5;
            end
            pA(ind) = pA(ind)+41;
        end
    end
    
    % scale the current levels to match the scaling of the measured levels
    
    xx = 0:400;
    [y1,x] = hist(pA,0:5:400);
    y1 = y1/sum(y1);
    y1 = min(y1,0.03);
    [y2,~] = hist(measured_levels,x);
    y2 = y2/sum(y2);
    y2 = min(y2,0.03);
    y1 = interp1(x,y1,0:400);
    y2 = interp1(x,y2,0:400);
    
    function value = fun(m,b,y1,y2)
        yy1 = interp1((0:400)*m+b,y1,0:400,'linear','extrap'); % scale the levels and reinterpolate
        value = sum((y2-yy1).^2);
    end
    
    params = fminsearch(@(p) fun(p(1),p(2),y1,y2), [range(measured_levels)/range(pA), min(measured_levels)-min(pA)*range(measured_levels)/range(pA)]);
    scale = params(1);
    offset = params(2);
    
    pA = pA * scale + offset;
    pA_std = pA_std * scale;

end