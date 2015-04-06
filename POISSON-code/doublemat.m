function [ A ] = doublemat( A )
    A = reshape(repmat(A',[2,1]),[2*numel(A) 1]);
end

