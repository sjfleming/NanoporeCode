function p = fit_levels(file, num, howtofit, pstart)
% fit_levels fits the Oxford levels to a physical model
% Stephen Fleming
% 2017/09/26
    
    %%    

    if nargin<4

        pstart = [0.3130
        0.7635
       -2.4443
       -0.9542
        2.3219
       -0.1784
       -0.1928
       -2.4485
        0.0832
        2.7365
       -0.1458
        0.2808
       -1.5316
       -0.7043
        2.1009
        0.3481
        1.1094
       -3.1951
       -0.6563
        2.3939
        3.1218
        0.5943
       55.3071
       70.1762]; % stretch: pA / base

    end

    %%
    
    d = load(file);
    list = arrayfun(@(x) d.model_data{num}.kmer(:,x)', 1:size(d.model_data{num}.kmer,2), 'uniformoutput', false); % list of 5-mers
    x = d.model_data{num}.level_mean;
    
    %%
    
    % fitting using matlab's built in optimization functionality
    
    if strcmp(howtofit,'matlab')

        op = optimset('maxfunevals', 31*100, 'maxiter', 31*100, 'display', 'iter');
        p = fminsearch(@(p) sum((x' - current_model_physical_2(list,p,0)).^2), pstart, op);

        p = normalizep(p);
    
    end
    
    %%
    
    if strcmp(howtofit,'sgd')
        
        rate = 0.01;
        iterations = 1000;
        for k = 1:iterations
            listing = randsample(numel(list),numel(list));
            for i = 1:numel(listing)
                kmer = list(listing(i));
                % difference from target
                delQ = (current_model_fun_sophisticated(kmer,p,0) - x(listing(i)));
                % update parameters
                pFa = p(1:5)' - mean(p(1:5));
                pFc = p(6:10)' - mean(p(6:10));
                pFg = p(11:15)' - mean(p(11:15));
                pFt = p(16:20)' - mean(p(16:20));
                seq = nt2int(list{i});
                Alogic = seq==1;
                Clogic = seq==2;
                Glogic = seq==3;
                Tlogic = seq==4;
                F_sum = sum(Alogic.*pFa + Clogic.*pFc + Glogic.*pFg + Tlogic.*pFt);
                % HERE
            end
            figure(2)
            plot(p','o')
            drawnow;
        end
        
        p = normalizep(p);
    
    end

    %%

    f = figure(1);
    clf
    plot(x,x,'ok')
    hold on
    plot(x,current_model_physical_2(list,pstart,0),'o')
    plot(x,current_model_physical_2(list,p,0),'o')
    dcm_obj = datacursormode(f);
    set(dcm_obj,'UpdateFcn',@(~,obj) list{obj.DataIndex})
    
    %%
    function p = normalizep(p)
        p(1:5) = p(1:5) - mean(p(1:5));
        p(6:10) = p(6:10) - mean(p(6:10));
        p(11:15) = p(11:15) - mean(p(11:15));
        p(16:20) = p(16:20) - mean(p(16:20));
    end
    
end