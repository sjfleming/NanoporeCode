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

%%

% plot bases in pore versus interaction position, at given voltage

v = 200; % voltage of interest
e = 5; % energy of interest
f = 18;
p = 9:17; % positions of interest

xx = [];
yy = cell(1,numel(p));

% look through files
for i = 1:numel(files)
    
    d = load(['data/' files{i}]);
    if ismember(d.position,p) && ...
            d.energy == e && ...
            d.force == f && ...
            d.voltage == v
        
        ind = find(d.position==xx,1,'first'); % locate x value
        if isempty(ind)
            xx(end+1) = d.position;
            ind = numel(xx);
        end
        
        % compile data
        granularity = 100;
        % interpolation, z as a function of base
        zz = linspace(1,d.mc.n,granularity);
        for j = 1:numel(d.mc.coordinates)
            z_locs = arrayfun(@(x) interp1((0:size(d.mc.coordinates{j},1)-1)*d.mc.l_k+d.mc.l_b,d.mc.coordinates{j}(:,3),d.mc.l_b*x,'spline'), zz);
            zz_ind_con = find(z_locs<=-7,1,'first');
            zz_ind_full = find(z_locs<=-8,1,'first');
            b(j) = zz(zz_ind_con); % base in pore constriction
            n(j) = zz(zz_ind_full); % bases within pore
        end
        
        % add y data to cell array (data are already thinned)
        yy{ind} = [yy{ind}, n];
        
    end
    
end


%%
figure(10)
hold on
errorbar(xx,cellfun(@(x) mean(x), yy),cellfun(@(x) std(x)/sqrt(numel(x)), yy),'ok')

% gaussian derivative
ft = fittype('a/sigma * (x-mu)/sigma^2 * exp(-(x-mu)^2/(2*sigma^2)) + b','coefficients',{'a', 'b', 'mu', 'sigma'},'independent',{'x'});
g = fit(xx',cellfun(@(x) mean(x), yy)',ft,'startpoint',[0.1,mean(yy{1}),13,2], ...
    'weights',cellfun(@(x) std(x)/sqrt(numel(x)), yy), 'lower', [0,0,0,0]);

% % cauchy derivative
% ft2 = fittype('a/sigma^2 * (x-mu) * (1+((x-mu)/sigma)^2)^-2 + b','coefficients',{'a', 'b', 'mu', 'sigma'},'independent',{'x'});
% g2 = fit(xx',cellfun(@(x) mean(x), yy)',ft2,'startpoint',[0.1,mean(yy{1}),13,2], ...
%     'weights',cellfun(@(x) std(x)/sqrt(numel(x)), yy), 'lower', [0,0,0,0.5]);

hold on
xxx = linspace(p(1),p(end),100);
plot(xxx,g(xxx),'--')
% plot(xxx,g2(xxx),'--')

ylabel('Number of nucleotides in pore')
xlabel('Nucleotide that interacts with constriction')
title(['Nucleotides in pore, ' num2str(f/100) 'pN/mV, ' num2str(e) 'kT interaction, ' num2str(v) 'mV'])
