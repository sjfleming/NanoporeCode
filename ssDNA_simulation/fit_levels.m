function p = fit_levels(file, num, modelfun, howtofit, pstart)
% fit_levels fits the Oxford levels to a physical model
% Stephen Fleming
% 2017/09/26
    
    %%    

    if nargin<5

        pstart = [-2.1141
    3.3616
   -3.7372
    1.5922
    0.8975
   -2.0784
    3.2771
   -2.3697
    1.6423
   -0.4712
   -0.3707
   -1.7028
    0.5570
    0.6200
    0.8966
    0.0828
   -0.8239
    3.5479
   -1.9386
   -0.8682
    0.5429
    0.3623
    0.2731
    0.2384
    0.5958
    0.5005
    0.2423
    0.2534
    0.4950
    0.7882
    0.4064
    0.5562
    0.1915
    0.3607
    0.6199
    0.4445
    0.8189
    0.5847
    0.6915
    0.6944
    0.4450
    0.0718
    0.2906];

% pstart = [0.2212
%     2.2573
%    -0.8391
%    -1.8038
%     0.1644
%    -1.0487
%    -0.4786
%    -0.9775
%     1.2076
%     1.2972
%    -0.8435
%     1.4279
%     0.9974
%    -1.1481
%    -0.4337
%    -2.2861
%     0.8477
%     2.2861
%    -0.1528
%    -0.6948
%     2.4998
%     1.4206
%    18.2357
%    25.3598];

    end

    %%
    
    d = load(file);
    list = arrayfun(@(x) d.model_data{num}.kmer(:,x)', 1:size(d.model_data{num}.kmer,2), 'uniformoutput', false); % list of 5-mers
    x = d.model_data{num}.level_mean;
    
    %%
    
    % fitting using matlab's built in optimization functionality
    
    if strcmp(howtofit,'matlab')

        op = optimset('maxfunevals', 31*100, 'maxiter', 31*100, 'display', 'iter');
        p = fminsearch(@(p) sum((x' - modelfun(list,p,0)).^2), pstart, op);

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
                delQ = (modelfun({kmer},p,0) - x(listing(i)));
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
            cost = sum((x' - modelfun(list,p,0)).^2);
            disp(num2str(cost))
        end
        
        p = normalizep(p);
    
    end
    
    %%
    
    figure(2)
    clf
    line([0 25],[0 0],'color','k')
    line([5.5 5.5],[-5 5],'color',[0.5,0.5,0.5])
    line(5+[5.5 5.5],[-5 5],'color',[0.5,0.5,0.5])
    line(10+[5.5 5.5],[-5 5],'color',[0.5,0.5,0.5])
    line(15+[5.5 5.5],[-5 5],'color',[0.5,0.5,0.5])
    
    if strcmp(howtofit,'sa')
        
        p = pstart;
        Tstart = 1000;
        Tend = 0.1;
        factor = 10;
        iterations = 100;
        
        currentcost = sum((modelfun(list,p,0) - x').^2);
        
        % annealing schedule
        T = logspace(log10(Tstart),log10(Tend),iterations);
        
        for k = 1:iterations
            
            inds = randsample(numel(p),numel(p));
            acc = 0;
            for i = 1:numel(inds)
                
                % proposal
                pr = p;
                thisind = inds(i);
                pairind = [];
                if thisind <=20  % only pair half the time
                    if rand()<0.5
                        % a random index from same group of five but not itself
                        pairind = (floor((thisind-1)/5))*5 + randsample(find(1:5~=(rem(thisind-1,5)+1)),1);
                    end
%                     continue; % don't mess with forces
                else
                    continue;
                end
                % pr(inds(i)) = pr(inds(i)) + randn(1)*pr(inds(i))/factor;
                movement = randn(1)/factor;%*max(0.1,pr(thisind));
                pr(thisind) = pr(thisind) + movement;
                pr(pairind) = pr(pairind) - movement;
                
                % cost
                cost = sum((modelfun(list,pr,0) - x').^2);
                
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
            disp(['T = ' num2str(T(k)) ', cost = ' num2str(currentcost) ', acceptance: ' num2str(round(acc/numel(p)*100)) '%'])
        end
        
        p = normalizep(p);
    
    end

    %%

    f = figure(1);
    clf
    plot(x,x,'ok')
    hold on
    plot(x,modelfun(list,pstart,0),'o')
    plot(x,modelfun(list,p,0),'o')
    dcm_obj = datacursormode(f);
    set(dcm_obj,'UpdateFcn',@(~,obj) list{obj.DataIndex})
    drawnow;
    
    %%
    function p = normalizep(p)
%         p(1:5) = p(1:5) - mean(p(1:5));
%         p(6:10) = p(6:10) - mean(p(6:10));
%         p(11:15) = p(11:15) - mean(p(11:15));
%         p(16:20) = p(16:20) - mean(p(16:20));
        p(21:end) = abs(p(21:end));
        % enforce homopolymer levels
        
%         % modify the total force
%         for ine = 1:4
%             indices = (1:5)+(ine-1)*5;
%             if abs(sum(p(indices))/20.48)/0.5 > 1 % dx in bases
%                 p(indices) = p(indices) / (abs(sum(p(indices))/20.48)/0.5);
%             end
%         end
        
%         % change the resistances
%         p(21:end) = abs(p(21:end));
%         p(21:25) = p(21:25) + (140/52.5778 - (sum(p(21:25)) + 5.921*exp(-((sum(p(1:5))/20.48)+140*0.18)/22.35) - 1.906))/5; % A 52.5778
%         p(26:30) = p(26:30) + (140/52.4166 - (sum(p(26:30)) + 5.921*exp(-((sum(p(6:10))/20.48)+140*0.18)/22.35) - 1.906))/5; % C 52.4166
%         p(31:35) = p(31:35) + (140/54.7827 - (sum(p(31:35)) + 5.921*exp(-((sum(p(11:15))/20.48)+140*0.18)/22.35) - 1.906))/5; % G 54.7827
%         p(36:40) = p(36:40) + (140/36.0750 - (sum(p(36:40)) + 5.921*exp(-((sum(p(16:20))/20.48)+140*0.18)/22.35) - 1.906))/5; % T 36.0750
%         p(21:end) = abs(p(21:end));
    end
    
end