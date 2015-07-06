function A_new = rshape( A, m, n )
%RSHAPE Smart replacement for Matlab's reshape function

    if nargin == 1
        % only given a single element, turn it into a column vector
        A_new = reshape(A, [numel(A) 1]);
    elseif nargin == 2
        % given A, m
        
        if numel(m) == 1
            % only given one dimension, calculate the second
            if mod(numel(A),m) ~= 0
                error('Dimension mismatch!');
            end
            m = [m,numel(A)/m];
        end
        
        A_new = reshape(A,m);
    else
        A_new = reshape(A,m,n);
    end
end

