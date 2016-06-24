function [errors,als] = find_errors(seq1,seq2)
% returns array of indices of errors between seq1 and seq2
% (indices into both arrays)

    [~,als] = seqalign(seq1,seq2);
    errors = als(any(~[1 1; diff(als)],2),:);

end