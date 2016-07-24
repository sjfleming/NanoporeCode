function [pA, pA_std, std_mean] = get_model_levels_oxford(seq, measured_levels, open_pore, voltage, temp)
% get_model_levels_oxford(seq, measured_levels) returns current levels for
% sequence 'seq' that are predicted by Oxford Nanopore Technologies, Inc.
% and scaled to mesh with the empirical CDF of the levels passed in,
% 'measured_levels'.
% Stephen Fleming, 2/8/16

    % Takes in a sequence in the form of 'ATCGTTCAAAGC...' with 'R'
    % denoting abasic sites

    % generate "states"
    states = get_states(nt2int(seq), 5); % k=5 for k-mer

    % get the current from Oxford's models
    load('models.mat')
    model = 1;
    pA = nan(size(states));
    std_mean = ones(size(states));
    % real parts
    pA(imag(states)==0) = model_data{model}.level_mean(states(imag(states)==0));
    pA_std(imag(states)==0) = model_data{model}.level_stdv(states(imag(states)==0));
    std_mean(imag(states)==0) = model_data{model}.sd_mean(states(imag(states)==0));
    
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
        end
    end
    
    % imaginary (abasic) parts
    %levs = model_data{model}.level_mean;
    %pA(imag(states)~=0) = imag(states(imag(states)~=0))*(max(levs)-min(levs))/2 + min(levs);
    
    % scale the current levels to match the scaling of the measured levels
    top = mean(measured_levels(measured_levels>open_pore*(0.44+(voltage-120)/1000)));
    if isnan(top)
        top = max(measured_levels);
        display('Trouble doing level scaling!!!!! (top level)')
    end
%     polyT = open_pore * ((voltage-120)*abs(temp-37)/5000*[1,1] + [0.08, 0.1]); % corrects for nonlinear IV in a heuristic way...
%     m = mode(round(measured_levels(measured_levels>polyT(1) & measured_levels<polyT(2)))); % trying to find poly T level
    sorted = sort(measured_levels,'ascend');
    m = mode(round(sorted(1:round(10/numel(seq)*numel(measured_levels)))));
    bottom = mean(measured_levels(measured_levels<m+2 & measured_levels>m-2));
    if isnan(bottom)
        bottom = open_pore*0.19;
        display('Trouble doing level scaling!!!!! (bottom level)')
    end
    %display(['[' num2str(bottom) ', ' num2str(top) ']'])
    pA = (pA-min(pA))*(top-bottom)/(range(pA)) + bottom;
    pA_std = pA_std*(top-bottom)/(range(pA));

end