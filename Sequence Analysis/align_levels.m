function [ dpath, best ] = align_levels( d0, d1, sd, si )
%ALIGN_LEVELS Align levels in two datasets
%   Uses Dynamic Time Warping (DTW) to align levels

    if (nargin < 3)
        % 0.4 basically means 2 standard devs. to replace a match
        % with two insertions. oughta keep at this level, i reckon
        sd = 2.0;
        si = 0.4;
    end

    % How many elements?
    n0 = numel(d0);
    n1 = numel(d1);

    scores = zeros(n0+1,n1+1);
    steps = zeros(n0+1,n1+1);

    scores(1,1) = 0;
    steps(1,1) = 3;

    % score returns something 0...1
    scr = @(i,j,sd) exp(-0.5*abs(d0(i)-d1(j)).^2/(sd)^2);

    % 1 is from above, 2 is from the left, 3 is diag
    scores(2:end,1) = (1:n0)*si;
    scores(1,2:end) = (1:n1)'*si;
    steps(2:end,1) = 1;
    steps(1,2:end) = 2;

    for i=1:n0
        for j=1:n1
            % check up and to the left etc
            sm = scr(i,j,sd);
            ss = [scores(i,j+1)+si, scores(i+1,j)+si, scores(i,j)+sm];
            
            [smax,ind] = max(ss);
            scores(i+1,j+1) = smax;
            steps(i+1,j+1) = ind;
        end
        if mod(i,100) == 0
            disp(num2str(i))
        end
    end

    % Now reconstruct the sequence that led here, as indices

    dpath = zeros(n0+n1,2);

    i = n0+1;
    j = n1+1;
    k = 1;
    while i>0 && j>0
        dpath(k,1) = i;
        dpath(k,2) = j;

        st = steps(i,j);
        switch (st)
            case 1
                i = i - 1;
            case 2
                j = j - 1;
            case 3
                i = i - 1;
                j = j - 1;
        end
        k = k + 1;
    end
    k = k - 1;

    fprintf('Aligned with %d total steps\n',k);
    dpath = flipud(dpath(1:k,:));
    
    best = scores(end,end);
    
end

