function pA = get_model_levels_oxford(seq, measured_levels, open_pore)
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
    % real parts
    pA(imag(states)==0) = model_data{model}.level_mean(states(imag(states)==0));
    
    % for now, input the levels exactly from measured data (2016_01_21_0009 at 1185sec has all these levels)
    if sum(imag(states)~=0)>0 % if there are imaginary numbers, we have abasics at those spots
        % go through each spot and match it to a known current level
        for ind = nonzeros((1:numel(states)).*(imag(states)~=0))'
            kmer = seq(ind:ind+4);
            switch kmer
                case 'RRRRR'
                    pA(ind) = 105.7;
                case 'RRRRT'
                    pA(ind) = 105.0;
                case 'RRRTT'
                    pA(ind) = 94.4;
                case 'RRTTT'
                    pA(ind) = 56.2;
                case 'RTTTT'
                    pA(ind) = 37.8;
                case 'GCGAR'
                    pA(ind) = 59.3;
                case 'CGARR'
                    pA(ind) = 72.7;
                case 'GARRR'
                    pA(ind) = 86.5;
                case 'ARRRA'
                    pA(ind) = 93.0;
                case 'RRRAC'
                    pA(ind) = 81.0;
                case 'RRACA'
                    pA(ind) = 62.0;
                case 'RACAT'
                    pA(ind) = 49.6;
                otherwise
                    display('Unknown abasic k-mer in sequence: cannot predict current level!');
            end
        end
    end
    
    % imaginary (abasic) parts
    %levs = model_data{model}.level_mean;
    %pA(imag(states)~=0) = imag(states(imag(states)~=0))*(max(levs)-min(levs))/2 + min(levs);
    
    % scale the current levels to match the scaling of the measured levels
    sorted_levels = sort(measured_levels,'descend');
    top = mean(measured_levels(measured_levels>open_pore*0.5));
    if isnan(top)
        top = sorted_levels(1);
    end
    m = mode(round(sorted_levels(end-10:end)));
    bottom = mean(measured_levels(measured_levels<m+1));
    %display(['[' num2str(bottom) ', ' num2str(top) ']'])
    pA = (pA-min(pA))*(top-bottom)/(range(pA)) + bottom;

end