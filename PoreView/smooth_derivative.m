function [deriv] = smooth_derivative( d, n )
%SMOOTH_DERIVATIVE is a filter for the data which finds the derivative
%   but in a smooth way.  Something...
%   Stephen Fleming, June 23, 2015
    
    %% Find differences between the mean for n points after and n points before each point
    
    after = zeros(n,numel(d));
    before = zeros(n,numel(d));
    
    for i = 1:n
        after(i,:) = [d(i+1:end)' zeros(1,i)];
        before(i,:) = [zeros(1,i) d(1:end-i)'];
    end
    
    deriv = mean(after,1) - mean(before,1);
    deriv(1:n) = zeros(1,n);
    deriv((end-n+1):end) = zeros(1,n);
    
    deriv = deriv';
    
end