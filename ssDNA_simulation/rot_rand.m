function R = rot_rand(d)
    % Get randomized rotation matrix, using Arvo 1992

    th = d*(rand()-0.5)*2*pi;
    phi = rand()*2*pi;
    z = d*rand()*2;
    r = sqrt(z);
    Vx = sin(phi)*r;
    Vy = cos(phi)*r;
    Vz = sqrt(2-z);
    st = sin(th);
    ct = cos(th);
    Sx = Vx*ct-Vy*st;
    Sy = Vx*st+Vy*ct;
    
    R = [-Vx*Sx+ct, -Vx*Sy+st, Vx*Vz;...
             -Vy*Sx-st, -Vy*Sy+ct, Vy*Vz;...
             -Vz*Sx,    -Vz*Sy,    1-z];

end
