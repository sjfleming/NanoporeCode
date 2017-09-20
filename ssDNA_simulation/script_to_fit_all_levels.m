%%

% find(strcmp('GGGGG',arrayfun(@(x) model_data{1}.kmer(:,x), 1:size(model_data{1}.kmer,2), 'uniformoutput', false)),1,'first')

%%

model_data{1}.level_mean([1 342 683 1024])
% A 52.5778
% C 52.4166
% G 54.7827
% T 36.0750

%%

list = arrayfun(@(x) model_data{1}.kmer(:,x), 1:size(model_data{1}.kmer,2), 'uniformoutput', false);

x = model_data{1}.level_mean;
y = zeros(size(x));

% p = fminsearch(@(p) sum((x' - currents(list,p)).^2), 0.2*ones(1,5));
% disp(p/sum(p))

op = optimset('maxfunevals', 15*10000, 'maxiter', 15*10000);
p = fminsearch(@(p) sum(abs(x' - current_model_fun(list,p))), [0.2*ones(1,5), ...
    -1, -2, 1, 2, 1,    -10, -1, 0, 1, 10], op);

figure(1)
clf
plot(x,x,'ok')
hold on
plot(x,current_model_fun(list,[ones(1,5), zeros(1,10)]),'o')
plot(x,current_model_fun(list,p),'o')

%%

seq = 'TTTTTGGGAAATTTTTAAAAACCCCCACCCCC';
inds = arrayfun(@(y) find(strcmp(seq(y:y+4),arrayfun(@(x) model_data{1}.kmer(:,x), 1:size(model_data{1}.kmer,2), 'uniformoutput', false)),1,'first'), 1:(numel(seq)-4));

figure(2)
plot(1:(numel(seq)-4),model_data{1}.level_mean(inds),'sk-')
hold on
plot(1:(numel(seq)-4),current_model_fun(arrayfun(@(x) model_data{1}.kmer(:,x), inds, 'uniformoutput', false), [p(1:5), zeros(1,10)]),'o--')
plot(1:(numel(seq)-4),current_model_fun(arrayfun(@(x) model_data{1}.kmer(:,x), inds, 'uniformoutput', false), p),'o--')

set(gca,'xtick',1:(numel(seq)-4),'xticklabelrotation',90,'xticklabel',arrayfun(@(x) seq(x:x+4), 1:(numel(seq)-4), 'uniformoutput',false))

%%

list = arrayfun(@(x) model_data{1}.kmer(:,x), 1:size(model_data{1}.kmer,2), 'uniformoutput', false);

x = model_data{1}.level_mean;
y = zeros(size(x));

op = optimset('maxfunevals', 67*100, 'maxiter', 67*100, 'display', 'iter');

pstart = ...
    [0.1305
    0.7407
    0.9199
    0.1385
         0.01
    1.0305
    3.9122
    5.6369
    1.5993
    2.7900
    4.9751
    6.2979
   -1.5192
   -0.0473
    2.6470
    4.1913
   -2.3632
   -1.1266
    0.6979
    1.6433
   -0.3373
    1.1492
    4.0669
    5.6421
    0.9656
    2.2866
    4.5371
    5.9129
   -3.8149
   -2.4656
    0.4166
    2.1354
   -4.6802
   -3.2674
   -1.3527
   -0.0051
    0.6450
    1.8498
    2.7776
   -3.9245
   -2.6207
   -0.8871
    0.0983
   -4.1092
   -3.0107
   -1.3185
   -0.1606
   -6.3345
   -4.7292
   -3.9113
   -2.8460
   -1.6483
   -0.6721
    0.7578
    1.3790
   -4.0356
   -2.9878
   -1.0518
   -0.2017
   -2.5383
   -1.2927
    0.7761
    1.7983
   -4.1019
   -3.7850
   -0.8033
   -0.0450];

p = fminsearch(@(p) sum((x' - current_model_fun_sophisticated(list,p)).^2), pstart, op);
% p = fminsearch(@(p) sum((x' - current_model_fun_sophisticated(list,p)).^2), SGD_init, op);
% p = fminsearch(@(p) sum((x' - current_model_fun_sophisticated(list,p)).^2), ...
%     [0.1214;    0.3613;    0.3920;    0.1136;    0.017; (p_A-p_not_A)-mean(p_A-p_not_A); (p_G-p_not_G)-mean(p_G-p_not_G)], op);

f = figure(2);
clf
plot(x,x,'ok')
hold on
plot(x,current_model_fun_sophisticated(list,[ones(1,5), zeros(1,64)]),'o')
plot(x,current_model_fun_sophisticated(list,p),'o')

dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn',@(~,obj) list{obj.DataIndex}')

%%

[num2str(p,'%10.4f'), ['w';'w';'w';'w';'w'; regexprep(regexprep(string(dec2bin(1:31)),'1','A'),'0','*'); regexprep(regexprep(string(dec2bin(1:31)),'1','G'),'0','*')]]

%%

seq = 'TTTTTGGGAAATTTTTAAAAACCCCCACCCCC';
inds = arrayfun(@(y) find(strcmp(seq(y:y+4),arrayfun(@(x) model_data{1}.kmer(:,x), 1:size(model_data{1}.kmer,2), 'uniformoutput', false)),1,'first'), 1:(numel(seq)-4));

figure(3);
clf
plot(1:(numel(seq)-4),model_data{1}.level_mean(inds),'sk-')
hold on
plot(1:(numel(seq)-4),current_model_fun_sophisticated(arrayfun(@(x) model_data{1}.kmer(:,x), inds, 'uniformoutput', false), [p(1:5); zeros(64,1)]),'o--')
plot(1:(numel(seq)-4),current_model_fun_sophisticated(arrayfun(@(x) model_data{1}.kmer(:,x), inds, 'uniformoutput', false), SGD_params),'o--')

set(gca,'xtick',1:(numel(seq)-4),'xticklabelrotation',90,'xticklabel',arrayfun(@(x) seq(x:x+4), 1:(numel(seq)-4), 'uniformoutput',false))


%%

% param estimation from data

p_G = cellfun(@(z) mean(arrayfun(@(a) model_data{1}.level_mean(a), z)), ...
    arrayfun(@(x) find(strcmp(cellfun(@(y) regexprep(y','[^G]','*'), list, 'uniformoutput',false),x)), ...
    regexprep(regexprep(string(dec2bin(1:31)),'1','G'),'0','*'), 'uniformoutput',false));

p_not_G = mean(arrayfun(@(z) model_data{1}.level_mean(z), ...
    find(strcmp(cellfun(@(y) regexprep(y','[^G]','*'), list, 'uniformoutput',false),'*****'))));

p_A = cellfun(@(z) mean(arrayfun(@(a) model_data{1}.level_mean(a), z)), ...
    arrayfun(@(x) find(strcmp(cellfun(@(y) regexprep(y','[^A]','*'), list, 'uniformoutput',false),x)), ...
    regexprep(regexprep(string(dec2bin(1:31)),'1','A'),'0','*'), 'uniformoutput',false));

p_not_A = mean(arrayfun(@(z) model_data{1}.level_mean(z), ...
    find(strcmp(cellfun(@(y) regexprep(y','[^A]','*'), list, 'uniformoutput',false),'*****'))));

%%

SGD_init = ...
    [0.1305
    0.7407
    0.9199
    0.1385
    0.0000
    0.9686
    4.2282
    5.7505
    1.5502
    2.7848
    5.5322
    6.7369
   -1.4061
    0.0092
    3.1960
    5.1370
   -2.4336
   -0.7890
    0.6733
    0.7251
   -0.4911
    0.8399
    4.3833
    5.7430
    0.7460
    1.8983
    5.3538
    6.0335
   -3.3116
   -2.0795
    0.2481
    2.5884
   -4.0524
   -2.8505
   -1.2935
   -0.0070
    0.7092
    1.8566
    3.0195
   -3.9402
   -2.8479
   -0.7622
    0.3417
   -4.3493
   -2.7644
   -0.4587
   -0.1206
   -6.0368
   -5.3601
   -3.5641
   -3.3628
   -1.8888
   -0.7159
    0.7329
    1.3812
   -4.1715
   -2.6234
   -0.8520
   -0.0507
   -2.3206
   -1.1721
    1.4108
    1.7287
   -4.0530
   -4.1268
   -0.7553
   -0.0543];

%%
% stochastic gradient descent

iterations = 1000;
rate = 0.5 / numel(list);
SGD_params = SGD_init;% [0.15; 0.3; 0.3; 0.1; 0.05; p_A - mean(p_A); p_G - mean(p_G)];
Q = @(p) sum((x' - current_model_fun_sophisticated(list,p)).^2);

allcombosA = regexprep(regexprep(string(dec2bin(1:31)),'1','A'),'0','*');
allcombosG = regexprep(regexprep(string(dec2bin(1:31)),'1','G'),'0','*');

figure(5)
clf
hold on
set(gca,'xtick',1:67,'xticklabelrotation',90,'xticklabel',{'w1','w2','w3','w4','w5',allcombosA{:},allcombosG{:}})
line([1 67],[0 0],'color','k')

for k = 1:iterations
    
    listing = randsample(numel(list),numel(list));
    for i = 1:numel(listing) % each example in "training set" (k-mer)
        
        % which parameters make a difference for this k-mer?
        param_ind_A = [];
        param_ind_G = [];
        
        % the corrections for As
        strA = regexprep(list{listing(i)}','[^A]','*'); % keeps only the As, rest are *
        if ~strcmp(strA,'*****')
            param_ind_A = find(strcmpi(strA,allcombosA),1,'first');
        end

        % the corrections for Gs
        strG = regexprep(list{listing(i)}','[^G]','*'); % keeps only the Gs, rest are *
        if ~strcmp(strG,'*****')
            param_ind_G = find(strcmpi(strG,allcombosG),1,'first');
        end
        
        % update params
        delQ = (current_model_fun_sophisticated(list(listing(i)),SGD_params) - model_data{1}.level_mean(listing(i)));
        if ~isempty(param_ind_A)
            SGD_params(5 + param_ind_A) = SGD_params(5 + param_ind_A) - rate * delQ;
        end
        if ~isempty(param_ind_G)
            SGD_params(5 + 31 + param_ind_G) = SGD_params(5 + 31 + param_ind_G) - rate * delQ;
        end
        
%         % update weights
%         if delQ<0 % we undershot the data
%             acg_inds = regexp(list{listing(i)}','[ACG]');
%             if ~isempty(acg_inds)
%                 rand_acg_ind = acg_inds(randi(numel(acg_inds))); % pick an A C or G at random
%                 vals = model_data{1}.level_mean([1 342 683 1024]);
%                 switch list{listing(i)}(rand_acg_ind)
%                     case 'A'
%                         val = vals(1);
%                     case 'C'
%                         val = vals(2);
%                     case 'G'
%                         val = vals(3);
%                     case 'T'
%                         val = vals(4);
%                 end
%                 SGD_params(rand_acg_ind) = SGD_params(rand_acg_ind) - rate * delQ * 0.01 * abs(model_data{1}.level_mean(listing(i)) - val)/val;
%             end
%         else
%             t_inds = regexp(list{listing(i)}','[T]');
%             if ~isempty(t_inds)
%                 rand_t_ind = t_inds(randi(numel(t_inds))); % pick an A C or G at random
%                 vals = model_data{1}.level_mean([1 342 683 1024]);
%                 switch list{listing(i)}(rand_acg_ind)
%                     case 'A'
%                         val = vals(1);
%                     case 'C'
%                         val = vals(2);
%                     case 'G'
%                         val = vals(3);
%                     case 'T'
%                         val = vals(4);
%                 end
%                 SGD_params(rand_t_ind) = SGD_params(rand_t_ind) - rate * delQ * 0.01 * abs(model_data{1}.level_mean(listing(i)) - val)/val;
%             end
%         end
%         SGD_params(1:5) = SGD_params(1:5)/sum(SGD_params(1:5));
        
    end
    
    disp([num2str(k) ', delQ = ' num2str(delQ) ', Q = ' num2str(Q(SGD_params))])
    figure(5)
    plot(1:numel(SGD_params),SGD_params,'o')
    drawnow
end

%%

g = figure(4);
clf
plot(x,x,'ok')
hold on
plot(x,current_model_fun_sophisticated(list,[ones(5,1); zeros(62,1)]),'o')
plot(x,current_model_fun_sophisticated(list,[SGD_params(1:5); zeros(62,1)]),'o')
plot(x,current_model_fun_sophisticated(list,SGD_params),'o')

dcm_obj = datacursormode(g);
set(dcm_obj,'UpdateFcn',@(~,obj) list{obj.DataIndex}')
