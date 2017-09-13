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

voltages = 60:20:200;
positions = 18:19;%9:17;
n = cell(numel(positions),numel(positions));
b = cell(numel(positions),1);
thinning = 2500;
independent_samples = 1000;
energy = 5; % factor of kT that interaction contributes
force = 18; % pN / 100mV
total = [];

for j = 1:numel(voltages)
    
    for f = 1:numel(positions)

        disp(['p = ' num2str(positions(f)) ', v = ' num2str(voltages(j))])
        
        for k = 1:round(independent_samples/(1e5/thinning))
            
            % run the sampling routine
            mc = ssDNA_MCMC('bases',28,'fixed_points',{1,[0,0,0]}, ...
                'force_function',@(d) force*(voltages(j)/100)*d(3), ...
                'boundary',@np_bnd,'initial_coordinates',init, ...
                'interaction_function',@(d) m2_constriction_interaction(d,positions(f),0.5,4.1*5,0.5,1.5));
            mc.run(1e5);
            
            % thin
            if isempty(total)
                mc.thin(thinning);
                total = mc;
            else
                mc.thin(thinning);
                total.coordinates = [total.coordinates, mc.coordinates];
            end
            
        end
        
        % save thinned data
        counter = 0;
        savename = ['data/mcmc_' sprintf('%04d',counter) '.mat'];
        while exist(savename,'file')~=0
            counter = counter+1;
            savename = ['data/mcmc_' sprintf('%04d',counter) '.mat'];
        end
        position = positions(f);
        voltage = voltages(j);
        mc = total;
        save(savename,'mc','position','voltage','energy','force','init','thinning');
        total = [];
        
    end
    
end