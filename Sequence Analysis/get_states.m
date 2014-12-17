function seqstates = get_states(seq, k)
% Get states from a sequence, using given kmer length (assumed 5)

    if (nargin < 2)
        k = 5;
    end
    
    seqstates = conv(double(seq)-1,4.^(0:k-1),'valid')+1;    
end