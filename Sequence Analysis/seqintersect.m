function [over1, over2, match] = seqintersect(seq1, seq2)

    [~, p, st] = swalign(seq1, seq2, 'Alphabet', 'NT');
    
    seqn1 = sum(p(1,:) ~= '-');
    seqn2 = sum(p(3,:) ~= '-');
    
    over1 = seq1(st(1):st(1)+seqn1-1);
    over2 = seq2(st(2):st(2)+seqn2-1);
    
    match = sum(p(2,:) == '|')/size(p,2);

end