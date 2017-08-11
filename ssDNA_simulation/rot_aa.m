function R = rot_aa(axis, th)
    % Rotation matrix about an axis by an angle
    x = axis(1);
    y = axis(2);
    z = axis(3);
    C = cos(th);
    S = sin(th);
    R = [x^2+(1-x^2)*C, (1-C)*x*y-z*S, (1-C)*x*z+y*S;...
        (1-C)*x*y+z*S, y^2+(1-y^2)*C, (1-C)*y*z-x*S;...
        (1-C)*x*z-y*S, (1-C)*y*z+x*S, z^2+(1-z^2)*C];
    
end
