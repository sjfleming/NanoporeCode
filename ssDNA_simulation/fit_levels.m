function p = fit_levels(file, num, pstart)
% fit_levels fits the Oxford levels to a physical model
% Stephen Fleming
% 2017/09/26
    
% %%    
%     pstart = ...
%     [0.0673 %weights
%     0.3819
%     0.4743
%     0.0714
%     0.0052
%     -0.5 % A
%     -1
%     -0.75
%     -0.2
%     -0.1
%     5 % A scale
%     0.1 % C
%     0.1
%     0
%     -0.1
%     -0.1
%     2 % C scale
%     -1 % G
%     -0.9
%     -0.5
%     0.2
%     0.3
%     4 % G scale
%     0.1 % T
%     0.2
%     0.2
%     0.2
%     0.1
%     1 % T scale
%     27 % stretch: pA / base
%     31]; % compress
% 
%     plower = [zeros(5,1); -1*ones(5,1); 0; -1*ones(5,1); 0; -1*ones(5,1); 0; -1*ones(5,1); 0; 0; 0];
%     pupper = [ones(5,1); ones(5,1); 10; ones(5,1); 10; ones(5,1); 10; ones(5,1); 10; 100; 100];

if nargin<3

% pstart = [0.0673 %weights
%     0.3819
%     0.4743
%     0.0714
%     0.0052
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
    
    %logic = cellfun(@(x) all(arrayfun(@(y) strcmp(y,'A') || strcmp(y,'G') || strcmp(y,'C'), x)), list);
    %list = list(logic);
    
    x = d.model_data{1}.level_mean;
    y = zeros(size(x));

    op = optimset('maxfunevals', 31*100, 'maxiter', 31*100, 'display', 'iter');
    p = fminsearch(@(p) sum((x' - current_model_physical_2(list,p,0)).^2), pstart, op);
    
    p(1:5) = p(1:5) - mean(p(1:5));
    p(6:10) = p(6:10) - mean(p(6:10));
    p(11:15) = p(11:15) - mean(p(11:15));
    p(16:20) = p(16:20) - mean(p(16:20));

    %%

    f = figure(1);
    clf
    plot(x,x,'ok')
    hold on
    plot(x,current_model_physical_2(list,pstart,0),'o')
    plot(x,current_model_physical_2(list,p,0),'o')
    dcm_obj = datacursormode(f);
    set(dcm_obj,'UpdateFcn',@(~,obj) list{obj.DataIndex})
    
end