%%

function script_plot_model_sequence(p)

% check how the model does on a sequence

% have a 'p' from a fit to current_model_physical_2

d = load('models.mat');
num = 1;
list = arrayfun(@(x) d.model_data{num}.kmer(:,x)', 1:size(d.model_data{num}.kmer,2), 'uniformoutput', false); % list of 5-mers
x = d.model_data{1}.level_mean;

seq = 'TTTTTGGGAAATTTTTAAAAACCCCCACCCCC';
inds = arrayfun(@(y) find(strcmp(seq(y:y+4),arrayfun(@(x) d.model_data{1}.kmer(:,x)', 1:size(d.model_data{1}.kmer,2), 'uniformoutput', false)),1,'first'), 1:(numel(seq)-4));

op = optimset('maxfunevals', 5*1000, 'maxiter', 5*1000);
fp = fminsearch(@(p) sum((x' - current_model_fun(list,p)).^2), [0.1305  0.7407  0.9199  0.1385  0.01], op);
fp = abs(fp)/sum(fp);

figure(3)
plot(1:(numel(seq)-4),d.model_data{1}.level_mean(inds),'sk-')
hold on
plot(1:(numel(seq)-4),current_model_fun(arrayfun(@(x) d.model_data{1}.kmer(:,x), inds, 'uniformoutput', false), fp),'o--')
plot(1:(numel(seq)-4),current_model_physical_3(arrayfun(@(x) d.model_data{1}.kmer(:,x)', inds, 'uniformoutput', false), p, 0),'o--')

set(gca,'xtick',1:(numel(seq)-4),'xticklabelrotation',90,'xticklabel',arrayfun(@(x) seq(x:x+4), 1:(numel(seq)-4), 'uniformoutput',false))

figure(5)
plot(x,x,'ok')
hold on
plot(x,current_model_fun(arrayfun(@(x) d.model_data{num}.kmer(:,x), 1:size(d.model_data{num}.kmer,2), 'uniformoutput', false),fp),'o')
plot(x,current_model_physical_3(list,p,0),'o')

end
