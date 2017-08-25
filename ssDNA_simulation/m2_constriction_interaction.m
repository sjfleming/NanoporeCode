function U = m2_constriction_interaction(coords, bases, length, energy, l_b, l_k)
% crude m2-mspa nanopore constriction interaction energy function
% coords = [x1,y1,z1; x2,y2,z2; ...] each bead
% bases = [b1, ...] each base that has this interaction
% length is the distance from constriction over which this force acts
%       suggested value = l_b = 0.5nm, from a quick guess based on
%       Bhattacharya "Water mediates..." ACS Nano 2016, figures 4c and 1b
% energy is the energy of this interaction in pN*nm, >0
%       for example 2*kT is about 8.2 pN*nm at room temp
% contributes an interaction energy, or not, depending on base positioning
% Stephen Fleming
% 8/17/17
    
    U = 0; % initially nothing
    scaling = 2/length^2; % energy * scaling has units of force/length
    
    % constriction position
    zcon = -7;
    
    % position of bases in z, roughly
    z_locs = arrayfun(@(x) interp1((0:size(coords,1)-1)*l_k+l_b,coords(:,3),l_b*bases(x),'spline'), 1:numel(bases)); % this interpolates in z
    % NOTE: based on fooling around, i suspect an off-by-one here...
    
    % are they in the interaction zone?
    % potential looks like a cut-off parabola that dips below zero near
    % zcon
    for i = 1:numel(z_locs)
        if z_locs(i) < zcon+length && z_locs(i) > zcon-length
            % potential looks like a cut-off parabola that dips below zero near zcon
            U = U - 1/2 * energy * scaling * (length^2 - (z_locs(i) - zcon)^2);
            % U = U - force*length; % just a bonus energy for being in that range
        end
    end

end