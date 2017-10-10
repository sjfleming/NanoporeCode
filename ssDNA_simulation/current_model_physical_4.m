function y = current_model_physical_4(list, p, debug)
% currents = current_model_physical_4(list, p) takes inputs of list, a
% cell array of level kmers, and p, a set of parameters, and returns
% currents for each kmer
% Gaussians: three
% Stephen Fleming
% 2017/10/09

    % resistances for A
    pRa = abs(p(1:5)');
    
    % resistances for C
    pRc = abs(p(6:10)');
    
    % resistances for G
    pRg = abs(p(11:15)');
    
    % resistances for T
    pRt = abs(p(16:20)');
    
    % params for Gaussians
    g_a = p(21:25);
    g_c = p(26:30);
    g_g = p(31:35);
    g_t = p(36:40);
    m_a = p(41:45);
    m_c = p(46:50);
    m_g = p(51:55);
    m_t = p(56:60);
    
    % interpolate the resistances based on movement of the strand
    xq = 0:0.01:6;
    pRa_int = interp1(-4:1:10,[zeros(1,5) pRa zeros(1,5)],xq,'pchip');
    pRc_int = interp1(-4:1:10,[zeros(1,5) pRc zeros(1,5)],xq,'pchip');
    pRg_int = interp1(-4:1:10,[zeros(1,5) pRg zeros(1,5)],xq,'pchip');
    pRt_int = interp1(-4:1:10,[zeros(1,5) pRt zeros(1,5)],xq,'pchip');

    % initialize output
    y = zeros(size(list));
    stretch = zeros(size(list));
    
    % loop through each kmer queried
    for i = 1:numel(list)
        
        seq = nt2int(list{i});
        Alogic = seq==1;
        Clogic = seq==2;
        Glogic = seq==3;
        Tlogic = seq==4;
        
        % calculate quasipotential
        xx = -1:0.01:1;
        T = 25 + 273.15;
        kT = T * 1.38e-23  / (1e-21); % energy scale kT
        U_total = 1/2 * 39.7 * (xx*0.5).^2 / kT; % initialize to Ustretch
        g = Alogic.*g_a' + Clogic.*g_c' + Glogic.*g_c' + Tlogic.*g_t';
        m = Alogic.*m_a' + Clogic.*m_c' + Glogic.*m_c' + Tlogic.*m_t';
        for j = 1:numel(seq)
            U_total = U_total + g(j) * (normpdf(xx,m(j),0.5) - normpdf(0,0,0.5));
        end
            
        % find energy minimum and assume force is adequate to move it there
        [~,min_ind] = min(U_total);
        dx = xx(min_ind);
        k = 20.48; % effective spring constant from MCMC in pN/nm
        F_sum = k*dx;
        f = F_sum+140*0.18; % grand total at 140mV
        R = 5.921*exp(-f/22.35) - 1.906; % change from baseline 140mV resistance (1.906G) based on fitting blocked pore resistance versus force
        
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