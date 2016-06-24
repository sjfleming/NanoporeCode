function [subset,inds] = randsubset(arr,num)

    num = min(num, numel(arr));
    inds = randperm(numel(arr),num);
    subset = arr(inds);

end