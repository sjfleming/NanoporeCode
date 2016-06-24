function [pA, pA_std] = get_model_levels_laszlo_M2(seq, measured_levels)
% get_model_levels_laszlo_M2(seq, measured_levels) returns current levels for
% sequence 'seq' that are predicted by M2-MspA measurements by Laszlo et al
% and scaled to mesh with the empirical CDF of the levels passed in,
% 'measured_levels'.
% Stephen Fleming, 6/1/16

    % Takes in a sequence in the form of 'ATCGTTCAAAGC...' with 'R'
    % denoting abasic sites

    % generate "states"
    k = 4;
    states = get_states(nt2int(seq), k); % k=4 for k-mer

    % get the current from Oxford's models
    d = load('qmerdatabase.mat');
%     load('models.mat');
%     model = 1;
    pA = nan(size(states));
    % real parts
    pA(imag(states)==0) = arrayfun(@(x) x.mean, d.qmerdatabase(states(imag(states)==0))) * 133; % convert into something like pA
%     pA(imag(states)==0) = model_data{model}.level_mean(states);
    pA_std(imag(states)==0) = arrayfun(@(x) x.error, d.qmerdatabase(states(imag(states)==0))) * 133;
%     pA_std(imag(states)==0) = model_data{model}.level_stdv(states);
    
    % for now, input the levels exactly from measured data (2016_01_21_0009 at 1185sec has all these levels)
    if sum(imag(states)~=0)>0 % if there are imaginary numbers, we have abasics at those spots
        % go through each spot and match it to a known current level
        for ind = nonzeros((1:numel(states)).*(imag(states)~=0))'
            kmer = seq(ind:ind+k-1);
            switch kmer
                case 'RRRR'
                    pA(ind) = 69.2+15;
                    pA_std(ind) = 0.9;
                case 'RRRT'
                    pA(ind) = 61.6+15;
                    pA_std(ind) = 1.5;
                case 'RRTT'
                    pA(ind) = 36.3+5;
                    pA_std(ind) = 3.5;
                case 'RTTT'
                    pA(ind) = 27.4+5;
                    pA_std(ind) = 1.1;
                case 'CGAR'
                    pA(ind) = 39.2;
                    pA_std(ind) = 1.5;
                case 'GARR'
                    pA(ind) = 41.6;
                    pA_std(ind) = 1.7;
                case 'ARRR'
                    pA(ind) = 57.1;
                    pA_std(ind) = 1.2;
                case 'RRRA'
                    pA(ind) = 63.9;
                    pA_std(ind) = 1.0;
                case 'RRAC'
                    pA(ind) = 56.3;
                    pA_std(ind) = 1.1;
                case 'RACA'
                    pA(ind) = 47.1;
                    pA_std(ind) = 1.5;
                otherwise
                    display(['Unknown abasic k-mer in sequence: guessing!   ' kmer]);
                    pA(ind) = (sum(kmer=='R')*105 + sum(kmer=='A')*53 + sum(kmer=='C')*52 + sum(kmer=='G')*55 + sum(kmer=='T')*36)/10;
                    pA_std(ind) = 1.5;
            end
        end
    end
    
    % scale the current levels to match the scaling of the measured levels
    
    [cdf1,lev1] = ecdf(pA);
    [cdf2,lev2] = ecdf(measured_levels);
    
    function value = fun(m,b,cdf1,cdf2,lev1,lev2)
        value = sum((interp1(cdf2,lev2,linspace(0.05,0.95,10),'linear','extrap')'-(m.*interp1(cdf1,lev1,linspace(0.05,0.95,10),'linear','extrap')'+b)).^2);
    end
    
    params = fminsearch(@(p) fun(p(1),p(2),cdf1,cdf2,lev1,lev2), [std(lev1)/std(lev2), mean(lev1)-mean(lev2)]);
    scale = params(1);
    offset = params(2);
    
    pA = pA * scale + offset;
    pA_std = pA_std * scale;

end