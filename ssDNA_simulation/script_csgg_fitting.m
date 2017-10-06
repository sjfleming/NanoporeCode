%%

clear all
d = load('CsgG_model.mat');
kmers = arrayfun(@(x) d.CsgG_model.kmer(:,x)', 1:size(d.CsgG_model.kmer,2), 'uniformoutput', false);
current = d.CsgG_model.level_mean;

% cut down by the ones with no weight

logic = (d.CsgG_model.weight ~= 0) & (current > 40);
current = current(logic);
current_sd = d.CsgG_model.sd_mean(logic);
kmers = kmers(logic);

% the re-factored and normalized data

x = cell2mat(cellfun(@(x) [nt2int(x)==1, nt2int(x)==2, nt2int(x)==3, nt2int(x)==4], kmers, 'uniformoutput', false)');
% x = x - mean(x(:,1));
% x = x / std(x(:,1));
% y = current - mean(current);
% y = y / std(y);
y = current;

%%

%clear all
d = load('models.mat');
kmers = arrayfun(@(x) d.model_data{1}.kmer(:,x)', 1:size(d.model_data{1}.kmer,2), 'uniformoutput', false); % list of 5-mers
current = d.model_data{1}.level_mean;

% the re-factored and normalized data

x2 = cell2mat(cellfun(@(x) [nt2int(x)==1, false, nt2int(x)==2, false, nt2int(x)==3, false, nt2int(x)==4, false], kmers, 'uniformoutput', false)');
y2 = current;
x = [x; x2];
y = [y; y2];

%%

nodes = 1;
w = randn(size(x,2),nodes)*0.5; % random initial weights
w = abs(w + 140/mean(current)/6); % fixing to be about right
s = 1.1;
b = 0.1;
cutpt = size(x,1)-size(x2,1)+1; % first index of MspA data

output = 140./(x*w);

epoch = 100;
eta = linspace(0.01,0.0001,epoch);

figure(1)
clf
plot(y,output,'o')

for i = 1:epoch
    ind = randsample(1:size(x,1),size(x,1));
    for j = 1:numel(ind)
        if ind(j)<cutpt
            dEdw = (y(ind(j))-output(ind(j))) * w / y(ind(j)) .* x(ind(j),:)';
            dEds = 0;
            dEdb = 0;
        else
            dEdw = (y(ind(j))-output(ind(j))) * (w*s+b) / y(ind(j)) .* x(ind(j),:)';
            dEds = sum((y(ind(j))-output(ind(j))) * w / y(ind(j)) .* x(ind(j),:)') / sum(x(ind(j),:));
            dEdb = (y(ind(j))-output(ind(j))) / (y(ind(j)) * sum(x(ind(j),:)));
        end
        w = w - eta(i)*dEdw;
        s = s - eta(i)*dEds;
        b = b - eta(i)*dEdb;
        w = abs(w);
        s = abs(s);
        b = abs(b);
        output = 140./([x(1:cutpt-1,:)*w; x(cutpt:end,:)*(w*s+b)]);
        %deriv = 140./w;
    end
    clf
    line([30 120],[30 120],'color','k')
    hold on
    plot(y,output,'o')
    ylim([30 110])
    xlim([30 110])
    title(num2str(i))
    drawnow
    disp(num2str(sum((y-output).^2)/numel(y)))
end

figure(2)
plot(w,'o')
line(size(x,2)/4*ones(1,2)+0.5,[0 1])
line(2*size(x,2)/4*ones(1,2)+0.5,[0 1])
line(3*size(x,2)/4*ones(1,2)+0.5,[0 1])

