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
    
    figure(2)
    clf
    
    if strcmp(howtofit,'sgd')
        
        p = pstart;
        rate = 0.0001;
        iterations = 1000;
        for k = 1:iterations
            listing = randsample(numel(list),numel(list));
            for i = 1:numel(listing)
                kmer = list{listing(i)};
                % difference from target
                delQ = (current_model_physical_2({kmer},p,0) - x(listing(i)));
                if delQ<=1
                    continue;
                end
                % update parameters
                pF(1,:) = p(1:5)' - mean(p(1:5));
                pF(2,:) = p(6:10)' - mean(p(6:10));
                pF(3,:) = p(11:15)' - mean(p(11:15));
                pF(4,:) = p(16:20)' - mean(p(16:20));
                seq = nt2int(kmer);
                Alogic = seq==1;
                Clogic = seq==2;
                Glogic = seq==3;
                Tlogic = seq==4;
                F_sum = sum(Alogic.*pF(1,:) + Clogic.*pF(2,:) + Glogic.*pF(3,:) + Tlogic.*pF(4,:));
                % HERE
                % pick one of the bases involved
                if numel(find([sum(Alogic)>0,sum(Clogic)>0,sum(Glogic)>0,sum(Tlogic)>0]))==1
                    base = find([sum(Alogic)>0,sum(Clogic)>0,sum(Glogic)>0,sum(Tlogic)>0]);
                else
                    base = randsample(find([sum(Alogic)>0,sum(Clogic)>0,sum(Glogic)>0,sum(Tlogic)>0]),1);
                end
                F_sum_base = sum((seq==base).*pF(base,:)); % sum of forces due to that random base
                % tweak either the forces or the current factor
                % when delQ<0, we want the current to go up
                % by stretching
                if sum(seq==base)==1
                    ind = (base-1)*5+find(seq==base);
                else
                    ind = (base-1)*5+randsample(find(seq==base),1);
                end
                if F_sum_base>0
                    p(ind) = p(ind) + delQ*39.7/2/p(end) * rate;
                else
                    p(ind) = p(ind) + delQ*39.7/2/p(end-1) * rate;
                end
                p((base-1)*5+1:base*5) = p((base-1)*5+1:base*5) - mean(p((base-1)*5+1:base*5)); % norm
                % do something about factors....
                if isnan(sum(p))
                    disp('')
                end
                %p = normalizep(p);
            end
            
            figure(2)
            hold on
            plot(p','o')
            drawnow;
%             disp(num2str(numel(p)))
            cost = sum((x' - current_model_physical_2(list,p,0)).^2);
            disp(num2str(cost))
        end
        
        p = normalizep(p);
    
    end
    
    %%
    
    figure(2)
    clf
    
    if strcmp(howtofit,'sa')
        
        p = pstart;
        Tstart = 10;
        iterations = 1000;
        
        currentcost = sum((current_model_physical_2(list,p,0) - x').^2);
        
        for k = 1:iterations
            
            % annealing schedule
            T = linspace(Tstart,1,iterations);
            
            inds = randsample(numel(p),numel(p));
            acc = 0;
            for i = 1:numel(inds)
                
                % proposal
                pr = p;
                pr(inds(i)) = pr(inds(i)) + randn(1)*pr(inds(i))/10;
                
                % cost
                cost = sum((current_model_physical_2(list,pr,0) - x').^2);
                
                % prob
                prob = exp(-(cost-currentcost)/T(k));
                if rand() < prob
                    % accept proposal and step
                    p = pr;
                    currentcost = cost;
                    acc = acc+1;
                end
                
                if isnan(sum(p))
                    disp('')
                end
                p = normalizep(p);
            end
            
            figure(2)
            hold on
            plot(p','o')
            drawnow;
            disp(['T = ' num2str(T(k)) ', cost = ' num2str(cost) ', acceptance: ' num2str(round(acc/numel(p)*100)) '%'])
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