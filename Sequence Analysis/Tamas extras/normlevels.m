function levels = normlevels(levels, levels1)
% normalize levels to 50 pA mean, 5 pA deviation, or to second set of
% levels (returning normalized first one)
    mval = 50;
    sval = 5;
    
    if (nargin == 2)
        mval = mean(levels1(:,1));
        sval = std(levels1(:,1));
    end

    %fprintf('scale: %f\n',sval/std(levels));
    
    levels = sval*levels/std(levels);
    %fprintf('shift: %f\n',mean(levels)-mval);
    levels = levels - mean(levels) + mval;

end