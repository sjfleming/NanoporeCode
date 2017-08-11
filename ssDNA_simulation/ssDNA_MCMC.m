classdef ssDNA_MCMC < handle
    % ssDNA_MCMC is a Markov chain Monte Carlo sampler of ssDNA position
    % configurations, inspired by Tamas Szalay's PoreMC.
    % Stephen Fleming
    % 8/10/17
    
    % freely jointed chain
    % U_{i,s} = 1/2 * k_s (l_{i+1} - l_i)^2 % stretching
    % U_{i,b} = -k_b * cos(theta_i) % bending
    
    properties
        
        % Kuhn length of ssDNA = 1.5nm
        l_k = 1.5;

        % bending energy of ssDNA
        k_b = 4; % pN*nm
        
        % stretching of ssDNA
        k_s = 800; % pN/nm, "Overstretching B-DNA...", Science 1996, -Smith, Cui, Bustamante, parameter 'S' on p. 798
        
        % number of bases of ssDNA
        n = 30;
        
        % length-per-base of ssDNA
        l_b = 0.5; % nm
        
        % contour length
        L = n*l_b; % nm
        
        % number of joints in the freely-jointed chain
        N = ceil(L/l_k); % next integer above contour length / Kuhn length
        
        % temperature
        T = 273.15 + 25; % Kelvin
        kT = T * 1.38e-23 / (4.11e-21); % k_B T, thermal energy, in pN*nm (=4.11e-21 Joules)
        
        % constraints
        fixed_points = cell(0,N);
        
        % all inputs from user
        in = [];
        
        % data
        coordinates = cell(0); % elements of cell array are time, each contains [x_1, y_1, z_1; x_2, y_2, z2; ...], in nm
        current_coords = []; % [x_1, y_1, z_1; x_2, y_2, z2; ...], in nm
        
    end
    
    methods (Static, Access = public)
    
        % constructor
        function obj = ssDNA_MCMC(varargin)
            % parse inputs
            obj.in = obj.parseInputs(varargin{:});
            obj.l_k = obj.in.l_k;
            obj.k_b = obj.in.k_b;
            obj.k_s = obj.in.k_s;
            obj.n = obj.in.n;
            obj.l_b = obj.in.l_b;
            obj.L = obj.n * obj.l_b;
            obj.N = ceil(obj.L/obj.l_k);
            obj.T = obj.in.T;
            obj.kT = obj.T * 1.38e-23  / (4.11e-21); % in pN*nm
            obj.fixed_points = obj.in.fixed_points;
            
            % initialize configuration of ssDNA
            if isempty(obj.coordinates)
                
            end
        end
        
        function run
            % run the mcmc sampler using metropolis-hastings
            
            % generate uniform probability r = U(0,1)
            r = rand;
            
            % generate proposed move
            test_coords = obj.propose();
            
            % caluclate energy of proposed move
            U = obj.energy(test_coords);
            
            % accept or reject, and take sample
            if r < U
                obj.accept();
            end
            obj.sample();
            
        end
        
        function test_coords = propose(obj)
            % propose a move
            
            
        end
        
        function U = energy(obj, coords)
            % calculate the energy of a proposed move
            U_s = 1/2 * obj.k_s * sum(sum(diff(coords,1).^2)); % stretching
            cosTheta = obj.calculateAngles(coords);
            U_b = -1 * obj.k_b * cosTheta;
            deltaU = U_s + U_b + obj.constraintEnergy(coords);
            U = exp(-deltaU / obj.kT);
        end
        
        function cosTheta = calculateAngles(obj,coords)
            % calculate angles of each segment with respect to previous one
            vectors = diff(coords,1); % now [dx1, dy1, dz1; dx2, dy2, dz2; ...]
            lengths = sum(vectors.^2,2); % vector lengths as [l1; l2; l3; ...]
            cosTheta = [0; vectors(1:end-1,:).*vectors(2:end,:)./lengths(1:end-1)./lengths(2:end); 0]; % N elements, first and last are zero
        end
        
        function U_c = constraintEnergy(obj, coords)
            % calculate the energy having to do with the constraints
            % imposed as boundary conditions
            
        end
        
        function accept(obj, test_coords)
            % accept proposal
            % move and update internal parameters
            obj.current_coords = test_coords;
        end
        
        function sample(obj)
            % take a sample
            obj.coordinates(:,end+1) = obj.current_coords;
        end
        
        function plotSnapshot(obj)
            
        end
        
        function plotOverlay(obj)
            
        end

        function in = parseInputs(varargin)
            % parse all inputs so all methods can use them easily
            p = inputParser;
            
            % set up the inputs
            addOptional(p, 'n', 30, @(x) x>6); % number of bases
            addOptional(p, 'l_k', 1.5, @(x) x>0.2); % Kuhn length of ssDNA, in nm
            addOptional(p, 'k_b', 4, @(x) x>0); % bending energy, in pN*nm
            addOptional(p, 'k_s', 800, @(x) x>0); % stretching modulus, in pN/nm
            addOptional(p, 'l_b', 0.5, @(x) x>0 && x<1); % length per base of ssDNA, in nm
            addOptional(p, 'T', 273.15 + 25, @(x) x>273.15 && x<373.15); % temperature in Kelvin
            addOptional(p, 'fixed_points', cell(0), @(x) iscell(x)); % {base_number_1, [x1, y1, z1], base_number_2, [x2,y2,z2], ...}
            addOptional(p, 'boundary', [], @(x) x>273.15 && x<373.15); % temperature in Kelvin

            obj.fixed_points = obj.in.fixed_points;
            
            % parse
            parse(p,varargin{:});
            in = p.Results;
        end
    
    end
    
end