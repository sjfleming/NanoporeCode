classdef ssDNA_MCMC < handle
    % ssDNA_MCMC is a Markov chain Monte Carlo sampler of ssDNA position
    % configurations, inspired by Tamas Szalay's PoreMC.
    % Stephen Fleming
    % 8/10/17
    
    % freely jointed chain
    % U_{i,s} = 1/2 * k_s ((l_{i+1} - l_i) - l_k)^2 % stretching
    % U_{i,b} = -k_b * cos(theta_i) % bending
    
    properties
        % tunable step size parameter, ratio of energies of steps to kT
        step;
        % Kuhn length of ssDNA = 1.5nm
        l_k;
        % bending energy of ssDNA, maybe 4 pN*nm
        k_b;
        % stretching of ssDNA
        k_s; % 800 pN/nm, "Overstretching B-DNA...", Science 1996, -Smith, Cui, Bustamante, parameter 'S' on p. 798
        % number of bases of ssDNA
        n;
        % length-per-base of ssDNA, 0.5nm
        l_b;
        % contour length, nm
        L;
        % number of joints in the freely-jointed chain, next integer above contour length / Kuhn length
        N;
        % temperature
        T; % Kelvin
        kT; % k_B T, thermal energy, in pN*nm (=4.11e-21 Joules)
        % constraints
        fixed_points;
        boundary;
        force_function;
        force_values;
        % initial ssDNA coordinates
        initial_coordinates;
        % all inputs from user
        in;
        % counting accpetance ratios
        proposal = '';
        count = struct();
        % keep track of energy
        current_energy = 0;
        % data
        coordinates = cell(0); % elements of cell array are time, each contains [x_1, y_1, z_1; x_2, y_2, z2; ...], in nm
        current_coords = []; % [x_1, y_1, z_1; x_2, y_2, z2; ...], in nm
        
    end
    
    methods (Access = public)
    
        % constructor
        function obj = ssDNA_MCMC(varargin)
            % parse inputs
            obj.in = obj.parseInputs(varargin{:});
            obj.k_b = obj.in.k_b;
            obj.k_s = obj.in.k_s;
            obj.n = obj.in.bases;
            obj.l_b = obj.in.l_b;
            obj.L = obj.n * obj.in.l_b;
            obj.N = ceil(obj.L/obj.in.l_k) + 1;
            obj.l_k = obj.L / (obj.N-1); % use the calculated length for segments
            obj.T = obj.in.T;
            obj.kT = obj.in.T * 1.38e-23  / (1e-21); % in pN*nm
            obj.fixed_points = obj.in.fixed_points;
            obj.boundary = obj.in.boundary;
            obj.force_function = obj.in.force_function;
            obj.force_values = obj.in.force_values;
            if ~isempty(obj.in.initial_coordinates)
                obj.initial_coordinates = obj.in.initial_coordinates;
            else
                obj.initial_coordinates = [(0:obj.N-1)'*obj.l_k, zeros(obj.N,2)];
            end
            obj.current_coords = obj.initial_coordinates;
            obj.step = obj.in.step;
            
            % initialize counters
            obj.count.proposed = struct();
            obj.count.proposed.translations = 0;
            obj.count.proposed.rotations = 0;
            obj.count.proposed.crankshafts = 0;
            obj.count.accepted = struct();
            obj.count.accepted.translations = 0;
            obj.count.accepted.rotations = 0;
            obj.count.accepted.crankshafts = 0;
            
            % initialize configuration of ssDNA
            if isempty(obj.coordinates)
                
            end
        end
        
        function run(obj, samples)
            % run the mcmc sampler using metropolis-hastings
            
            for i = 1:samples
            
                % generate uniform probability r = U(0,1)
                r = rand();

                % generate proposed move
                test_coords = obj.propose();

                % caluclate energy of proposed move
                U = obj.energy(test_coords);

                % accept or reject proposal, and either way take a sample
                if r < exp(-(U-obj.current_energy)/obj.kT)
                    obj.accept(test_coords, U);
                end
                obj.sample();
            
            end
            
        end
        
        function test_coords = propose(obj)
            % propose a move
            r = randi(3); % random integer: 1, 2, or 3
            if r==1
                obj.proposal = 'translation';
                obj.count.proposed.translations = obj.count.proposed.translations + 1;
                test_coords = obj.propose_translation();
            elseif r==2
                obj.proposal = 'rotation';
                obj.count.proposed.rotations = obj.count.proposed.rotations + 1;
                test_coords = obj.propose_rotation();
            else
                obj.proposal = 'crankshaft';
                obj.count.proposed.crankshafts = obj.count.proposed.crankshafts + 1;
                test_coords = obj.propose_crankshaft();
            end
            % ensure fixed points are fixed
            % this is probably the worst way to do it... work on that
            if ~isempty(obj.fixed_points)
                for i = 1:numel(obj.fixed_points)/2
                    test_coords(obj.fixed_points(2*i-1)) = obj.fixed_points(2*i);
                end
            end
        end
        
        function test_coords = propose_translation(obj)
            % propose a "translation" move
            % which beads are involved
            inds = sort(randi(obj.N,1,2)); % two integers from 1 to N, in order
            delta = randn(1,3);
            unit = delta/sqrt(sum(delta.^2));
            delta_scaled = unit * sqrt(obj.step * 2 * obj.kT / obj.k_s * pi / 2) * rand(1)*obj.step; % pi is for avg cos(theta)^2, randn is for no hard max
            test_coords = obj.current_coords;
            test_coords(inds(1):inds(2),:) = test_coords(inds(1):inds(2),:) + repmat(delta_scaled,diff(inds)+1,1);
        end
        
        function test_coords = propose_rotation(obj)
            % propose a "rotation" move
            % which beads are involved
            i = randi(obj.N); % random integer from 1 to N
            dtheta = tanh(sqrt(2/obj.k_b) * obj.step^2); % empirical
            % the maximum dtheta is 1, and min is 0
            R = rot_rand(dtheta); % random small angle rotation matrix from Tamas' implementation of Arvo 1992
            fixed_pt = obj.current_coords(i,:); % this is the fixed bead
            j = 1;
            if rand()<0.5 % rotate either all beads before or all beads after fixed bead
                j = obj.N;
            end
            inds = sort([i,j]); % bead numbers in order, one being fixed, other being either end
            vectors = obj.current_coords(inds(1):inds(2),:) - repmat(fixed_pt,diff(inds)+1,1); % vectors to each bead from fixed point
            new_vectors = vectors*R; % rotate the vectors
            test_coords = obj.current_coords;
            test_coords(inds(1):inds(2),:) = repmat(fixed_pt,diff(inds)+1,1) + new_vectors;
        end
        
        function test_coords = propose_crankshaft(obj)
            % propose a "crankshaft" move% which beads are involved
            % which beads are involved
            inds = sort(randi(obj.N,1,2)); % two integers from 1 to N, in order
            axis = obj.current_coords(inds(2),:) - obj.current_coords(inds(1),:); % vector between the two beads
            theta = rand * 2*pi; % random number on [0, 2*pi]
            R = rot_aa(axis, theta); % rotation maxtrix for rotating around an axis by angle theta
            fixed_pt = obj.current_coords(inds(1),:);
            vectors = obj.current_coords(inds(1):inds(2),:) - repmat(fixed_pt,diff(inds)+1,1); % vectors to each bead from fixed point
            new_vectors = vectors*R; % rotate the vectors
            test_coords = obj.current_coords;
            test_coords(inds(1):inds(2),:) = repmat(fixed_pt,diff(inds)+1,1) + new_vectors;
        end
        
        function U = energy(obj, coords)
            % calculate the energy of a proposed move
            vectors = diff(coords,1); % now [dx1, dy1, dz1; dx2, dy2, dz2; ...]
            lengths = sqrt(sum(vectors.^2,2)); % vector lengths as [l1; l2; l3; ...]
            U_s = 1/2 * obj.k_s * sum((lengths - obj.l_k).^2); % stretching
            cosTheta = obj.calculateAngles(coords);
            % U_b = -1 * obj.k_b * (sum(cosTheta) - sum(obj.calculateAngles(obj.initial_coordinates)));
            U_b = -1 * obj.k_b * sum(cosTheta);
            % deltaU = U_s + U_b + obj.constraintEnergy(coords) - obj.current_energy;
            U = U_s + U_b + obj.constraintEnergy(coords);
            % U = exp(-deltaU / obj.kT);
            % U = U_b / obj.kT;
        end
        
        function cosTheta = calculateAngles(obj,coords)
            % calculate angles of each segment with respect to previous one
            vectors = diff(coords,1); % now [dx1, dy1, dz1; dx2, dy2, dz2; ...]
            lengths = sqrt(sum(vectors.^2,2)); % vector lengths as [l1; l2; l3; ...]
            cosTheta = [0; sum(vectors(1:end-1,:).*vectors(2:end,:),2)./lengths(1:end-1)./lengths(2:end); 0]; % N elements, first and last are zero
        end
        
        function U_f = constraintEnergy(obj, coords)
            % calculate the energy having to do with the constraints
            % imposed as boundary conditions
            U_f = 0;
            if ~isempty(obj.boundary)
                
            end
            if ~isempty(obj.force_function)
                
            elseif ~isempty(obj.force_values)
                
            end
        end
        
        function accept(obj, test_coords, test_energy)
            % accept proposal
            % move and update internal parameters
            obj.current_coords = test_coords;
            obj.current_energy = test_energy;
            switch obj.proposal
                case 'translation'
                    obj.count.accepted.translations = obj.count.accepted.translations + 1;
                case 'rotation'
                    obj.count.accepted.rotations = obj.count.accepted.rotations + 1;
                case 'crankshaft'
                    obj.count.accepted.crankshafts = obj.count.accepted.crankshafts + 1;
            end
        end
        
        function sample(obj)
            % take a sample
            obj.coordinates{end+1} = obj.current_coords;
        end
        
        function reset_all(obj)
            % clear all the samples from this simulation
            obj.coordinates = cell(0);
            obj.current_coords = obj.initial_coordinates;
            obj.current_energy = 0;
            % clear counters
            obj.count.proposed.translations = 0;
            obj.count.proposed.rotations = 0;
            obj.count.proposed.crankshafts = 0;
            obj.count.accepted.translations = 0;
            obj.count.accepted.rotations = 0;
            obj.count.accepted.crankshafts = 0;
        end
        
        function overall = acceptance_ratios(obj)
            % display acceptance ratios
            overall = sum(structfun(@(x) x, obj.count.accepted)) ...
                / sum(structfun(@(x) x, obj.count.proposed))*100;
            disp(['Overall: ' num2str(overall) '% accepted'])
            disp(['Translations: ' num2str(obj.count.accepted.translations  ...
                / obj.count.proposed.translations * 100) '% accepted'])
            disp(['Rotations: ' num2str(obj.count.accepted.rotations  ...
                / obj.count.proposed.rotations * 100) '% accepted'])
            disp(['Crankshafts: ' num2str(obj.count.accepted.crankshafts  ...
                / obj.count.proposed.crankshafts * 100) '% accepted'])
        end
        
        function plot_snapshot(obj, time_index, fig)
            % plot a 3d line plot of ssDNA for given time index in figure
            figure(fig)
            hold on
            plot3(obj.coordinates{time_index}(:,1), ...
                obj.coordinates{time_index}(:,2), ...
                obj.coordinates{time_index}(:,3),'o-')
        end
        
        function plot_overlay(obj, thinning, fig)
            % plot a 3d line plot of ssDNA at (all) timepoints in figure
            % if thinning is 5, plots every fifth timepoint
            for i = 1:thinning:numel(obj.coordinates)
                obj.plot_snapshot(i, fig);
            end
        end

        function in = parseInputs(obj,varargin)
            % parse all inputs so all methods can use them easily
            p = inputParser;
            
            % set up the inputs
            addOptional(p, 'l_k', 1.5, @(x) x>0.2); % Kuhn length of ssDNA, in nm
            addOptional(p, 'k_b', 4, @(x) x>0); % bending energy, in pN*nm
            addOptional(p, 'k_s', 800, @(x) x>0); % stretching modulus, in pN/nm
            addOptional(p, 'l_b', 0.5, @(x) x>0 && x<1); % length per base of ssDNA, in nm
            addOptional(p, 'T', 273.15 + 25, @(x) x>273.15 && x<373.15); % temperature in Kelvin
            addOptional(p, 'fixed_points', cell(0), @(x) iscell(x) && numel(x{2})==3); % {base_number_1, [x1, y1, z1], base_number_2, [x2,y2,z2], ...}
            addOptional(p, 'boundary', [], @(x) islogical(feval(x,[0,0,0]))); % logical function: true inside boundary, false outside
            addOptional(p, 'force_function', [], @(x) isnumeric(feval(x,[0,0,0]))); % numeric function that gives force (pN) as a function of location (3d nm)
            addOptional(p, 'force_values', [], @(x) isnumeric(feval(x,[0,0,0]))); % values from which to interpolate force (pN)
            addOptional(p, 'initial_coordinates', [], @(x) all(isnumeric(x)) && size(x,2)==3); % starting ssDNA coordinates, for Kuhn segments: [x1,y1,z1; ...; xN,yN,zN] (nm)
            addOptional(p, 'bases', 30, @(x) x>6); % number of bases
            addOptional(p, 'step', 1, @(x) x>0); % step size, in an energy ratio to kT
            
            % parse
            parse(p,varargin{:});
            in = p.Results;
        end
    
    end
    
end