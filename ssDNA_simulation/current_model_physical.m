function y = current_model_physical(list, p, debug)
% currents = current_model_physical(list, p) takes inputs of list, a
% cell array of level kmers, and p, a set of parameters, and returns
% currents for each kmer
% Stephen Fleming
% 2017/09/26
    
    % weights
%     w = p(1:5); %[p(1:4) max(0,1-sum(p(1:4)))]; % find fifth weight given 4
%     w = w/sum(w); % enforce normalization
    w = [0.0673 %weights
    0.3819
    0.4743
    0.0714
    0.0052];
    
    % params for A
    pUa = p(1:5) - mean(p(1:5));
    
    % params for C
    pUc = p(6:10) - mean(p(6:10));
    
%     % params for G
%     pUg = p(16:20);
%     
%     % params for T
%     pUt = p(21:25);
%     
%     % params for stretching and compression's effect on current
%     current_stretch_ratio = p(26);
%     current_compress_ratio = p(27);
    
    current_stretch_ratio = 30;
    %current_compress_ratio = p(17);

    % initialize output
    y = zeros(size(list));
    
    % compute the stretching force given the quasi-potentials
    xx = -1:0.002:1; % in bases
    T = 25 + 273.15;
    kT = T * 1.38e-23  / (1e-21); % energy scale kT
    % parabola, Ustretch = 1/2 * (39.7pN/nm) * x^2, from simulation of base in constriction versus voltage
    % based on 2.4 base movement
    % with a force 23.8 pN
    % division by kT puts this into units of kT
    
    % loop through each kmer queried
    for i = 1:numel(list)
        
        Utotal = 1/2 * 39.7 * (xx*0.5).^2 / kT; % initialize to Ustretch
        
        for j = 1:numel(list{i})
        
            switch list{i}(j)
                case 'A'
                    pU = pUa';
                case 'C'
                    pU = pUc';
                case 'G'
                    pU = pUg';
                case 'T'
                    pU = pUt';
            end

            % create an interpolation around this point
            U = interp1(-4:1:10, [zeros(1,4) [pU(1) pU pU(end)] zeros(1,4)], xx+j, 'pchip');
            dU = U - U(ceil(numel(xx)/2)); % subtract value at center to get difference relative to unstretched
            Utotal = Utotal + dU; % add this, scaled, to running total quasi-potential
            if debug
                figure(3)
                plot(xx,dU)
                hold on
            end
            clear U dU

        end

        % now we have the full quasi-potential, Utotal
        [~,ind] = min(Utotal);
        if isempty(ind)
            disp('what...')
        end
        dx = xx(ind); % stretched by this much, positive being compression of strand
%         if dx<0
%             dI = current_stretch_ratio * dx;
%         else
%             dI = current_compress_ratio * dx;
%         end
        dI = current_stretch_ratio * dx;
        
        if debug
            figure(2)
            plot(xx,Utotal)
            hold on
        end
        
        % re-interpolate the weights based on movement of the strand
        w_shift = interp1(-4:1:10,[zeros(1,5) w' zeros(1,5)],(1:5)+dx,'pchip');
        w_shift = abs(w_shift)/sum(abs(w_shift)); % re-normalize
%         w_shift = w;
        
        % do a weighted average over bases
        for j = 1:numel(list{i})
            
            switch list{i}(j)
                case 'A'
                    y(i) = y(i) + 52.5778 * w_shift(j);
                case 'C'
                    y(i) = y(i) + 52.4166 * w_shift(j);
                case 'G'
                    y(i) = y(i) + 54.7827 * w_shift(j);
                case 'T'
                    y(i) = y(i) + 36.0750 * w_shift(j);
            end
            
        end
        
        y(i) = y(i) + dI; % add current from stretching
        
    end
    
    if debug
        disp('debug')
    end
    
end