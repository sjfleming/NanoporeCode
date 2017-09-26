function y = current_model_physical_2(list, p, debug)
% currents = current_model_physical_2(list, p) takes inputs of list, a
% cell array of level kmers, and p, a set of parameters, and returns
% currents for each kmer
% Stephen Fleming
% 2017/09/26
    
    % weights
%     w = p(1:5); %[p(1:4) max(0,1-sum(p(1:4)))]; % find fifth weight given 4
%     w = w/sum(w); % enforce normalization
%     w = [0.0673 %weights
%     0.3819
%     0.4743
%     0.0714
%     0.0052];
%     weightfun = fit([1:5]',w,'gauss1','startpoint',[2,2.5,1]);
    weightfun = @(x) normpdf(x,p(end-3),p(end-2));
    
    % params for A
    pFa = p(1:5)' - mean(p(1:5));
    
    % params for C
    pFc = p(6:10)' - mean(p(6:10));
    
    % params for G
    pFg = p(11:15)' - mean(p(11:15));
    
    % params for T
    pFt = p(16:20)' - mean(p(16:20));
    
    % params for stretching and compression's effect on current
    current_stretch_ratio = p(end-1);
    current_compress_ratio = p(end);

    % initialize output
    y = zeros(size(list));
    
    % compute the stretching force
    
    % loop through each kmer queried
    for i = 1:numel(list)
        
        seq = nt2int(list{i});
        Alogic = seq==1;
        Clogic = seq==2;
        Glogic = seq==3;
        Tlogic = seq==4;
        
        %Utotal = 1/2 * 39.7 * (xx*0.5).^2 / kT; % initialize to Ustretch
        k = 39.7;
        F_sum = sum(Alogic.*pFa + Clogic.*pFc + Glogic.*pFg + Tlogic.*pFt);
        
        % now we have the full quasi-potential, Utotal
        dx = 2*F_sum/k;
        if dx<0
            dI = -1 * current_stretch_ratio * dx;
        else
            dI = -1 * current_compress_ratio * dx;
        end
        
        % re-interpolate the weights based on movement of the strand
%         w_shift = interp1(-4:1:10,[zeros(1,5) w' zeros(1,5)],(1:5)+dx,'pchip');
        w_shift = weightfun((1:5)+dx);
        w_shift = abs(w_shift)/sum(abs(w_shift)); % re-normalize
%         w_shift = w';
        
        % do a weighted average over bases
        y(i) = sum((Alogic*52.5778 + Clogic*52.4166 + Glogic*54.7827 + Tlogic*36.0750) .* w_shift);
        y(i) = y(i) + dI; % add current from stretching
        
    end
    
    if debug
        disp('debug')
    end
    
end