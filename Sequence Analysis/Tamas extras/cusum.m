function v = cusum(v)

    csum = 0;
    for i=1:numel(v)
        csum = csum + v(i);
        if csum < 0
            csum = 0;
        end
        v(i) = csum;
    end

end