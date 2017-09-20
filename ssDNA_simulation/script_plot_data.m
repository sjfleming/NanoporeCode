% Stephen Fleming
% 9/12/17

% open and plot saved data
% generated from script_gen_data.m

%%

clear all

% find files

files = what('data/');
files = files.mat;
files = files(cellfun(@(x) strcmp(x(1:5),'mcmc_'), files));
cmap = parula(9);

%%

c = cmap(1,:);

% plot bases in pore versus interaction position, at given voltage

v = 160; % voltage of interest
e = 5; % energy of interest
f = 18;
p = 9:19; % positions of interest

xx = [];
yy = cell(1,numel(p));

% look through files
for i = 1:numel(files)
    
    d = load(['data/' files{i}]);
    % position = d.position;
    position = max(d.pos);
    if ismember(position,p) && ...
            d.energy == e && ...
            d.force == f && ...
            d.voltage == v
        
        ind = find(position==xx,1,'first'); % locate x value
        if isempty(ind)
            xx(end+1) = position;
            ind = numel(xx);
        end
        
        % compile data
        granularity = 1000;
        % interpolation, z as a function of base
        zz = linspace(d.mc.l_b/2,d.mc.n*d.mc.l_b,granularity);
        for j = 1:numel(d.mc.coordinates)
            z_locs = interp1((0:size(d.mc.coordinates{j},1)-1)*d.mc.l_k+d.mc.l_b/2, d.mc.coordinates{j}(:,3), zz,'spline');
            zz_ind_con = find(z_locs<=-7,1,'first');
            zz_ind_full = find(z_locs<=-8,1,'first');
            b(j) = zz(zz_ind_con)/d.mc.l_b; % base in pore constriction
            n(j) = zz(zz_ind_full)/d.mc.l_b; % bases within pore
        end
        
        % add y data to cell array (data are already thinned)
        yy{ind} = [yy{ind}, n];
        
    end
    
end

%%
figure(2)
hold on
errorbar(xx,cellfun(@(x) mean(x), yy),cellfun(@(x) std(x)/sqrt(numel(x)), yy),'o','color',c)

% gaussian derivative
ft = fittype('a/sigma * (x-mu)/sigma^2 * exp(-(x-mu)^2/(2*sigma^2)) + b','coefficients',{'a', 'b', 'mu', 'sigma'},'independent',{'x'});
g = fit(xx',cellfun(@(x) mean(x), yy)',ft,'startpoint',[0.1,mean(yy{1}),14,2], ...
    'weights',cellfun(@(x) std(x)/sqrt(numel(x)), yy), 'lower', [0,0,0,0], 'upper', [0.6, 18, 19, 4]);

% % cauchy derivative
% ft2 = fittype('a/sigma^2 * (x-mu) * (1+((x-mu)/sigma)^2)^-2 + b','coefficients',{'a', 'b', 'mu', 'sigma'},'independent',{'x'});
% g2 = fit(xx',cellfun(@(x) mean(x), yy)',ft2,'startpoint',[0.1,mean(yy{1}),13,2], ...
%     'weights',cellfun(@(x) std(x)/sqrt(numel(x)), yy), 'lower', [0,0,0,0.5]);

hold on
xxx = linspace(p(1),p(end),100);
plot(xxx,g(xxx),'-','color',c)
% plot(xxx,g2(xxx),'--')

display(['v = ' num2str(v) 'mV'])
g

ylabel('Number of nucleotides in pore')
xlabel('Nucleotide that interacts with constriction')
title(['Nucleotides in pore, ' num2str(f/100) 'pN/mV, ' num2str(e) 'kT interaction, ' num2str(v) 'mV'])

%%

Ib = @(x,p) p(1) * (p(2) - p(3) * x) * v;
p = [0.0976   23.6141    1.2579]; % based on a fit to open and blocked pore data from 20160629

figure(3)
I = Ib(cellfun(@(x) mean(x), yy),p);
errorbar(xx,I,Ib(cellfun(@(x) mean(x) - std(x)/sqrt(numel(x)), yy),p)-I,Ib(cellfun(@(x) mean(x) + std(x)/sqrt(numel(x)), yy),p)-I,'o','color',c)

set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0],'xtick',[11:16],'xticklabel',{'CCCCC','CCCCA','CCCAA','CCAAA','CAAAA','AAAAA'})
ylabel('Current (pA)')
title(['Predicted current, ' num2str(f/100) 'pN/mV, ' num2str(e) 'kT interaction, ' num2str(v) 'mV'])
ylim([51 82])
xlim([10.1 16.9])
