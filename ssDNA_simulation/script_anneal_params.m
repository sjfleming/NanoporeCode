%%

p = [2; 1; 0; -1; -2;
    0; 0; 0; 0; 0;
    -4; 3; 2; 0.5; -1.5;
    0; 0; 0; 0; 0;
    0.5*ones(15,1);
    0.9*ones(5,1)];

%%

p = [2; 1; 0; -1; -2;
    0; 0; 0; 0; 0;
    -4; 3; 2; 0.5; -1.5;
    0; 0; 0; 0; 0;
    0.4751
    0.6140
    0.7769
    0.5605
    0.5005
    0.4603
    0.5157
    0.6233
    0.5532
    0.5510
    0.4941
    0.6696
    0.8438
    0.5846
    0.5150
    0.4766
    0.5564
    0.6709
    0.6543
    0.5627];

%%

for i = 1:1000
    % simulated anneal, then hone in to a local minimum, then repeat bump
    % out and repeat
    p = fit_levels('models.mat',1,@(x,p,i) current_model_physical_3(x,p,i),'sa',p);%+randn(size(p))*0.05);
    p = fit_levels('models.mat',1,@(x,p,i) current_model_physical_3(x,p,i),'matlab',p);
end


%%

% attempt to get the params from the data

d = load('models.mat');
kmers = arrayfun(@(x) d.model_data{1}.kmer(:,x)', 1:size(d.model_data{1}.kmer,2), 'uniformoutput', false); % list of 5-mers
currents = d.model_data{1}.level_mean;

logicn = cellfun(@(x) ~isempty(regexp(x,'[^A][ACGT][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[A][ACGT][ACGT][ACGT][ACGT]')), kmers);
pa(1) = mean(currents(logic)) - mean(currents(logicn));
pas(1) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][^A][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][A][ACGT][ACGT][ACGT]')), kmers);
pa(2) = mean(currents(logic)) - mean(currents(logicn));
pas(2) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][^A][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][A][ACGT][ACGT]')), kmers);
pa(3) = mean(currents(logic)) - mean(currents(logicn));
pas(3) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][^A][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][A][ACGT]')), kmers);
pa(4) = mean(currents(logic)) - mean(currents(logicn));
pas(4) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][^A]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][A]')), kmers);
pa(5) = mean(currents(logic)) - mean(currents(logicn));
pas(5) = std(currents(logic));

logicn = cellfun(@(x) ~isempty(regexp(x,'[^C][ACGT][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[C][ACGT][ACGT][ACGT][ACGT]')), kmers);
pc(1) = mean(currents(logic)) - mean(currents(logicn));
pcs(1) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][^C][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][C][ACGT][ACGT][ACGT]')), kmers);
pc(2) = mean(currents(logic)) - mean(currents(logicn));
pcs(2) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][^C][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][C][ACGT][ACGT]')), kmers);
pc(3) = mean(currents(logic)) - mean(currents(logicn));
pcs(3) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][^C][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][C][ACGT]')), kmers);
pc(4) = mean(currents(logic)) - mean(currents(logicn));
pcs(4) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][^C]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][C]')), kmers);
pc(5) = mean(currents(logic)) - mean(currents(logicn));
pcs(5) = std(currents(logic));

logicn = cellfun(@(x) ~isempty(regexp(x,'[^G][ACGT][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[G][ACGT][ACGT][ACGT][ACGT]')), kmers);
pg(1) = mean(currents(logic)) - mean(currents(logicn));
pgs(1) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][^G][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][G][ACGT][ACGT][ACGT]')), kmers);
pg(2) = mean(currents(logic)) - mean(currents(logicn));
pgs(2) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][^G][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][G][ACGT][ACGT]')), kmers);
pg(3) = mean(currents(logic)) - mean(currents(logicn));
pgs(3) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][^G][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][G][ACGT]')), kmers);
pg(4) = mean(currents(logic)) - mean(currents(logicn));
pgs(4) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][^G]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][G]')), kmers);
pg(5) = mean(currents(logic)) - mean(currents(logicn));
pgs(5) = std(currents(logic));

logicn = cellfun(@(x) ~isempty(regexp(x,'[^T][ACGT][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[T][ACGT][ACGT][ACGT][ACGT]')), kmers);
pt(1) = mean(currents(logic)) - mean(currents(logicn));
pts(1) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][^T][ACGT][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][T][ACGT][ACGT][ACGT]')), kmers);
pt(2) = mean(currents(logic)) - mean(currents(logicn));
pts(2) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][^T][ACGT][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][T][ACGT][ACGT]')), kmers);
pt(3) = mean(currents(logic)) - mean(currents(logicn));
pts(3) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][^T][ACGT]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][T][ACGT]')), kmers);
pt(4) = mean(currents(logic)) - mean(currents(logicn));
pts(4) = std(currents(logic));
logicn = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][^T]')), kmers);
logic = cellfun(@(x) ~isempty(regexp(x,'[ACGT][ACGT][ACGT][ACGT][T]')), kmers);
pt(5) = mean(currents(logic)) - mean(currents(logicn));
pts(5) = std(currents(logic));

p = [-pa'-mean(-pa);
    -pc'-mean(-pc);
    -pg'-mean(-pg);
    -pt'-mean(-pt);
    0.5*ones(5,1);
    0.5*ones(5,1);
    0.5*ones(5,1);
    0.9*ones(5,1)];

%%

p = [0.4751
    0.6140
    0.7769
    0.5605
    0.5005
    0.4603
    0.5157
    0.6233
    0.5532
    0.5510
    0.4941
    0.6696
    0.8438
    0.5846
    0.5150
    0.4766
    0.5564
    0.6709
    0.6543
    0.5627
    -0.1*ones(5,1)+randn(5,1)*0.05
    randn(5,1)*0.1
    randn(5,1)*0.1
    0.1*ones(5,1)+randn(5,1)*0.1
    randn(20,1)*0.1];

%%

for i = 1:1000
    % simulated anneal, then hone in to a local minimum, then repeat bump
    % out and repeat
    p = fit_levels('models.mat',1,@(x,p,i) current_model_physical_4(x,p,i),'sa',p);%+randn(size(p))*0.05);
    p = fit_levels('models.mat',1,@(x,p,i) current_model_physical_4(x,p,i),'matlab',p);
end
