function inside = np_bnd(coord)
% crude nanopore boundary function
% coord = [x,y,z]
% true inside boundary
% false outside
% Stephen Fleming
% 8/16/17
    
    % pore radius
    R = 0.6; % nm
    
    % z location of constriction
    zcon = -7; % nm
    
    inside = true;
    
    % nanopore only exists from z = -8 to z = 0
    if coord(3) < 0 && coord(3) > -8
        
        % z = -7 is the narrowest constriction

        % from z = -8 to z = -7
        if coord(3) < -7
            if coord(1)^2 + coord(2)^2 >= R^2 - (coord(3)-zcon)*3*R^2
                inside = false;
            end

        % from z = -7 to z= 0
        else
            if coord(1)^2 + coord(2)^2 >= R^2 + (coord(3)-zcon)*(15/7)*R^2
                inside = false;
            end
        end
    
    end
    
end