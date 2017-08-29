%%

% plot the nanopore as three rings

hold on

R = 0.6;

x = linspace(-R,R,100);
y1 = sqrt(R^2-x.^2);
y2 = -sqrt(R^2-x.^2);
plot3([x,x],[y1,y2],ones(1,numel(x)*2)*-7,'k')
plot3([x,x],[y1,y2],ones(1,numel(x)*2)*-7,'k.')

x = linspace(-4*R,4*R,100);
y1 = sqrt((4*R)^2-x.^2);
y2 = -sqrt((4*R)^2-x.^2);
plot3([x,x],[y1,y2],ones(1,numel(x)*2)*0,'k')
plot3([x,x],[y1,y2],ones(1,numel(x)*2)*0,'k.')

x = linspace(-2*R,2*R,100);
y1 = sqrt((2*R)^2-x.^2);
y2 = -sqrt((2*R)^2-x.^2);
plot3([x,x],[y1,y2],ones(1,numel(x)*2)*-8,'k')
plot3([x,x],[y1,y2],ones(1,numel(x)*2)*-8,'k.')

%%

init = [0         0         0
   -0.0154    0.0235   -1.4465
   -0.3619    0.2136   -3.1522
   -0.1791    0.3043   -4.8435
   -0.0844    0.3333   -6.4507
    0.2508    0.8718   -7.8463
    0.5506    0.8994   -9.3131
    0.3292    1.0806  -10.8514
    0.9230    0.9985  -12.3720
    0.5799    1.9483  -13.6168
    1.2485    1.9650  -14.9118];

%%

% count number of bases in pore, from z = 0 to z = -8
% at various voltages

voltages = [60:20:180];
n = cell(numel(voltages),1);
thinning = 2500;

for f = 1:numel(voltages)

disp(['V = ' num2str(voltages(f))])
mc = ssDNA_MCMC('bases',30,'fixed_points',{1,[0,0,0]},'force_function',@(d) 18*(voltages(f)/100)*d(3),'boundary',@np_bnd,'initial_coordinates',init);
mc.run(100000);

n{f} = [];

for i = thinning:thinning:numel(mc.coordinates)
    
    j = find(mc.coordinates{i}(:,3)<-8,1,'first'); % index of first bead past pore
    granularity = 100;
    z = linspace(mc.coordinates{i}(j-1,3),mc.coordinates{i}(j,3),granularity); % list
    k = find(z<-8,1,'first');
    n{f}(end+1) = mc.l_k/mc.l_b * ((j-2) + (k-1)/granularity); % n = number of bases, first bead is the endpoint
    
end

end

%%

% plotting that

errorbar(voltages, (mean(n{1}) - cellfun(@(x) mean(x), n))/mean(n{1}) + 0.11, cellfun(@(x) std(x), n)/mean(n{1})/sqrt(40),'ok')
xlim([40 180])
set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
ylabel('fractional change in number of bases in pore')
xlabel('voltage (mV)')
title('ssDNA MCMC simulation')
ylim([0.08 0.3])

%%

% different way, using poisson ratio extension versus area thing, vol same

errorbar(voltages, (mean(n{1}) - cellfun(@(x) mean(x), n))/mean(n{1}) + 0.11, cellfun(@(x) std(x), n)/mean(n{1})/sqrt(40),'ok')
xlim([40 180])
set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
ylabel('fractional change in number of bases in pore')
xlabel('voltage (mV)')
title('ssDNA MCMC simulation')
ylim([0.08 0.3])

%% 

%%

% count number of bases in pore, from z = 0 to z = -8
% at various voltages
% with interactions

voltages = 60:20:180;
n = cell(numel(voltages),1);
b = cell(numel(voltages),1);
thinning = 2500;

for f = 1:numel(voltages)

disp(['V = ' num2str(voltages(f))])
mc = ssDNA_MCMC('bases',28,'fixed_points',{1,[0,0,0]},'force_function',@(d) 18*(voltages(f)/100)*d(3), ...
    'boundary',@np_bnd,'initial_coordinates',init,'interaction_function',@(d) m2_constriction_interaction(d,14,0.5,4.1*5,0.5,1.5)); % 5kT for base 14
mc.run(100000);

n{f} = [];
b{f} = [];

for i = thinning:thinning:numel(mc.coordinates)

    granularity = 1000;
    % interpolation, z as a function of base
    zz = linspace(1,mc.n,granularity);
    z_locs = arrayfun(@(x) interp1((0:size(mc.coordinates{i},1)-1)*mc.l_k+mc.l_b,mc.coordinates{i}(:,3),mc.l_b*x,'spline'), zz);
    zz_ind_con = find(z_locs<=-7,1,'first');
    zz_ind_full = find(z_locs<=-8,1,'first');
    b{f}(end+1) = zz(zz_ind_con); % base in pore constriction
    n{f}(end+1) = zz(zz_ind_full); % bases within pore
    
end

end

%% plot base in pore constriction

%figure()
errorbar(voltages,cellfun(@(x) mean(x), b),cellfun(@(x) std(x)/sqrt(numel(x)), b),'o')
set(gca,'fontsize',12,'outerposition',[0.01,0.01,0.98,0.98],'looseinset',[0,0,0,0])
xlim([40 200])
ylabel('base in pore constriction')
xlabel('voltage (mV)')
title('ssDNA MCMC simulation')

%%

% count number of bases in pore, from z = 0 to z = -8
% at various voltages
% with interactions

voltage = 160;
positions = 9:17;
n = cell(numel(voltages),1);
b = cell(numel(voltages),1);
thinning = 2500;
independent_samples = 100;

for f = 1:numel(positions)
    
    n{f} = [];
    b{f} = [];
    
    for k = 1:round(independent_samples/(1e5/thinning))
    
        disp(['p = ' num2str(positions(f))])
        mc = ssDNA_MCMC('bases',28,'fixed_points',{1,[0,0,0]},'force_function',@(d) 18*(voltage/100)*d(3), ...
            'boundary',@np_bnd,'initial_coordinates',init,'interaction_function',@(d) m2_constriction_interaction(d,positions(f),0.5,4.1*5,0.5,1.5)); % 5kT for base 14
        mc.run(1e5);

        for i = thinning:thinning:numel(mc.coordinates)

            granularity = 1000;
            % interpolation, z as a function of base
            zz = linspace(1,mc.n,granularity);
            z_locs = arrayfun(@(x) interp1((0:size(mc.coordinates{i},1)-1)*mc.l_k+mc.l_b,mc.coordinates{i}(:,3),mc.l_b*x,'spline'), zz);
            zz_ind_con = find(z_locs<=-7,1,'first');
            zz_ind_full = find(z_locs<=-8,1,'first');
            b{f}(end+1) = zz(zz_ind_con); % base in pore constriction
            n{f}(end+1) = zz(zz_ind_full); % bases within pore

        end

    end

end