function y = current_model_physical_3(list, p, debug)
% currents = current_model_physical_3(list, p) takes inputs of list, a
% cell array of level kmers, and p, a set of parameters, and returns
% currents for each kmer
% Stephen Fleming
% 2017/09/27
    
    % params for A
    pFa = p(1:5)' - mean(p(1:5));
    
    % params for C
    pFc = p(6:10)' - mean(p(6:10));
    
    % params for G
    pFg = p(11:15)' - mean(p(11:15));
    
    % params for T
    pFt = p(16:20)' - mean(p(16:20));
    
    % resistances for A
    pRa = abs(p(21:25)');
    
    % resistances for C
    pRc = abs(p(26:30)');
    
    % resistances for G
    pRg = abs(p(31:35)');
    
    % resistances for T
    pRt = abs(p(36:40)');
    
    % params for stretching and compression's effect on current
    %current_stretch_resistance_change_vs_dx = p(end-1);
%     min_vestibule_resistance = abs(p(end-2));
%     dRdx = abs(p(end-1));
%     range_vestibule_resistance = abs(p(end));
    
    % interpolate the resistances based on movement of the strand
    xq = 0:0.01:6;
    pRa_int = interp1(-4:1:10,[zeros(1,5) pRa zeros(1,5)],xq,'pchip');
    pRc_int = interp1(-4:1:10,[zeros(1,5) pRc zeros(1,5)],xq,'pchip');
    pRg_int = interp1(-4:1:10,[zeros(1,5) pRg zeros(1,5)],xq,'pchip');
    pRt_int = interp1(-4:1:10,[zeros(1,5) pRt zeros(1,5)],xq,'pchip');

    % initialize output
    y = zeros(size(list));
    
    % loop through each kmer queried
    for i = 1:numel(list)
        
        seq = nt2int(list{i});
        Alogic = seq==1;
        Clogic = seq==2;
        Glogic = seq==3;
        Tlogic = seq==4;
        
        k = 20.48; % effective spring constant from MCMC in pN/nm
        F_sum = sum(Alogic.*pFa + Clogic.*pFc + Glogic.*pFg + Tlogic.*pFt);
        
        % now we have the full quasi-potential, Utotal
        dx = F_sum/k; % in nm
        %dx = 0; % TESTING
        R = 0;
        %R = R + min_vestibule_resistance + range_vestibule_resistance./(1+exp(-dx/dRdx));
        f = F_sum+140*0.18; % grand total at 140mV
        R = R + 5.921*exp(-f/22.35) - 1.906; % change from baseline 140mV resistance (1.906G) based on fitting blocked pore resistance versus force
%         if dx<0
%             %R = R + current_stretch_resistance_change_vs_dx * dx;
%             R = R + current_compress_resistance_change_vs_dx * dx;
%         else
%             R = R + current_compress_resistance_change_vs_dx * dx;
%         end
        
        % add resistance of each base
        index_for_xq = -1*round((dx/0.5)/0.01) + (1:5)/0.01 + 1; % if dx>0, compression, so sample points x-dx
        index_for_xq = min(max(1,index_for_xq),numel(xq)); % make sure inds are in range (if they're out, they become zero resistance)
        R = R + sum(Alogic.*pRa_int(index_for_xq) + Clogic.*pRc_int(index_for_xq) + ...
            Glogic.*pRg_int(index_for_xq) + Tlogic.*pRt_int(index_for_xq));
        
        % calculate current y
        y(i) = 140 / R; % in pA, when R is in Gohm, and current is at 140mV
        stretch(i) = dx;
    end
    
    if debug
        disp('debug')
    end
    
end