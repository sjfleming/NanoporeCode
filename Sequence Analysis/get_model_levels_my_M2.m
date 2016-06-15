function [pA, pA_std] = get_model_levels_my_M2(seq, measured_levels, varargin)
% get_model_levels_oxford(seq, measured_levels) returns current levels for
% sequence 'seq' that are predicted by Oxford Nanopore Technologies, Inc.
% and scaled to mesh with the empirical CDF of the levels passed in,
% 'measured_levels'.
% Stephen Fleming, 6/2/16

    % Takes in a sequence in the form of 'ATCGTTCAAAGC...' with 'R'
    % denoting abasic sites

    % generate "states"
    states = get_states(nt2int(seq), 5); % k=5 for k-mer

    % get the current from Oxford's models
    d = load('models.mat');
    model = 1;
    pA = nan(size(states));
    % real parts
    pA(imag(states)==0) = d.model_data{model}.level_mean(states(imag(states)==0));
    pA_std(imag(states)==0) = d.model_data{model}.level_stdv(states(imag(states)==0));
    
    % for now, input the levels exactly from measured data (2016_01_21_0009 at 1185sec has all these levels)
    if sum(imag(states)~=0)>0 % if there are imaginary numbers, we have abasics at those spots
        % go through each spot and match it to a known current level
        for ind = nonzeros((1:numel(states)).*(imag(states)~=0))'
            kmer = seq(ind:ind+4);
            switch kmer
                case 'RRRRR'
                    pA(ind) = 105.3;
                    pA_std(ind) = 0.74;
                case 'RRRRT'
                    pA(ind) = 104.1;
                    pA_std(ind) = 0.64;
                case 'RRRTT'
                    pA(ind) = 93.0;
                    pA_std(ind) = 1.35;
                case 'RRTTT'
                    pA(ind) = 52.4;
                    pA_std(ind) = 5.79;
                case 'RTTTT'
                    pA(ind) = 38.1;
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
    
    % replace known levels for the M2 pore
    for ind = 1:numel(states)
        kmer = seq(ind:ind+4);
        switch kmer
            case 'TTTTT'
                pA(ind) = 35.9;
                pA_std(ind) = 0.75;
            case 'TTTTG'
                pA(ind) = 38.0;
                pA_std(ind) = 0.75;
            case 'TTTGG'
                pA(ind) = 40.5;
                pA_std(ind) = 1;
            case 'TTGGG'
                pA(ind) = 48.2;
                pA_std(ind) = 1.2;
            case 'TGGGA'
                pA(ind) = 58.5;
                pA_std(ind) = 1.5;
            case 'GGGAA'
                pA(ind) = 65.0;
                pA_std(ind) = 1.75;
            case 'GGAAA'
                pA(ind) = 67.3;
                pA_std(ind) = 1.5;
            case 'GAAAT'
                pA(ind) = 63.0;
                pA_std(ind) = 1.1;
            case 'AAATT'
                pA(ind) = 46.9;
                pA_std(ind) = 1.2;
            case 'AATTT'
                pA(ind) = 39.7;
                pA_std(ind) = 0.9;
            case 'ATTTT'
                pA(ind) = 37.5;
                pA_std(ind) = 0.8;
            case 'TTTTC'
                pA(ind) = 35.9;
                pA_std(ind) = 0.8;
            case 'TTTCG'
                pA(ind) = 42.8;
                pA_std(ind) = 0.8;
            case 'TTCGA'
                pA(ind) = 54.1;
                pA_std(ind) = 0.8;
            case 'TCGAT'
                pA(ind) = 60.7;
                pA_std(ind) = 0.8;
            case 'CGATC'
                pA(ind) = 49.8;
                pA_std(ind) = 0.8;
            case 'GATCA'
                pA(ind) = 50.6;
                pA_std(ind) = 0.8;
            case 'ATCAC'
                pA(ind) = 53.9;
                pA_std(ind) = 0.8;
            case 'TCACT'
                pA(ind) = 53.6;
                pA_std(ind) = 0.8;
            case 'CACTG'
                pA(ind) = 46.9;
                pA_std(ind) = 0.8;
            case 'ACTGG'
                pA(ind) = 47.5;
                pA_std(ind) = 0.8;
            case 'CTGGA'
                pA(ind) = 56.4;
                pA_std(ind) = 0.8;
            case 'TGGAA'
                pA(ind) = 61.3;
                pA_std(ind) = 0.8;
%             case 'GGAAC'
%                 pA(ind) = ;
%                 pA_std(ind) = 0.8;
            otherwise
                % nothing
        end
    end
    
%     figure()
%     plot(pA,'o-')
    
    % scale the current levels to match the scaling of the measured levels
    
    function value = fun(m,b,cdf1,cdf2,lev1,lev2)
        value = sum((interp1(cdf2,lev2,linspace(0.1,0.9,10),'linear','extrap')'-(m.*interp1(cdf1,lev1,linspace(0.1,0.9,10),'linear','extrap')'+b)).^2);
    end
    
    if isempty(varargin) || numel(varargin{1})==0
    
        cut = 0.75;
        [cdf1,lev1] = ecdf(pA(pA<cut*max(pA)));
        if isempty(measured_levels(measured_levels<cut*max(measured_levels)))
            display('Worning: not enough levels to construct a CDF at all... guessing.')
            scale = max(measured_levels)/max(pA);
            offset = 0;
        else
            [cdf2,lev2] = ecdf(measured_levels(measured_levels<cut*max(measured_levels)));
            params = fminsearch(@(p) fun(p(1),p(2),cdf1,cdf2,lev1,lev2), [std(lev1)/std(lev2), mean(lev1)-mean(lev2)]);
            scale = params(1);
            offset = params(2);
        end
        
    else
        
        scale = (varargin{1}{1}(2)-varargin{1}{1}(1))/(max(pA)-min(pA));
        offset = varargin{1}{1}(1) - min(pA)*scale;
        
    end
    
    pA = pA * scale + offset;
    pA_std = pA_std * scale;

end