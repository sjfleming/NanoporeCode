function [seq, merged] = seqmerge(seq1, seq2, nover)
    % Merge two sequences based on short overlap. The first sequence
    % should be the longer of the two for best performance.
    
    % length of overlap to use
    if nargin < 3
        nover = 50;
    end
    % if seq1 is short
    nover = min(nover, numel(seq1)-1);
    
    % get this length from the long one (seq1) on each end
    sleft = seq1(1:nover);
    sright = seq1(end-nover:end);
    
    % now do swalign on both and get the starting indices/p-s
    [~,pleft,stleft] = swalign(sleft,seq2,'Alphabet','NT');
    [~,pright,stright] = swalign(sright,seq2,'Alphabet','NT');
    
    % calculate the end indices
    eleft = stleft + sum(pleft([1,3],:)~='-',2) - 1;
    eright = stright + sum(pright([1,3],:)~='-',2) - 1;
    
    % and how many matches each one has
    mleft = sum(pleft(2,:) == '|');
    mright = sum(pright(2,:) == '|');
    
    % check if we have good enough matching, meaning at least 80% match
    % and that at least 50% of nover was used....
    accleft = mleft/size(pleft,2);
    accright = mright/size(pright,2);
    if size(pleft,2) < 0.5*nover
        accleft = 0;
    end
    if size(pright,2) < 0.5*nover
        accright = 0;
    end
    if max(accleft,accright) < 0.75
        fprintf(2,'Sequences do not match well enough to merge\n');
        %error('Sequences do not match well enough to merge');
        merged = false;
        seq = seq1;
        return
    end
    
    % now we use whichever has more matches
    if accleft > accright
        % mleft means seq2 is to the left
        % so eleft should be close to end of seq2
        seq = [seq2(1:stleft(2)-1) seq1];
    else
        % mright means seq2 is to the right
        % so stright should be close to the start of seq2
        seq = [seq1 seq2(eright(2)+1:end)];
    end
    merged = true;
    
end