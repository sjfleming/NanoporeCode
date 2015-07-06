function swtest()

    seq1 = int2nt(randi(4,[1,4000]));
    seq2 = seq1;
    seq2(randperm(numel(seq2),300)) = int2nt(randi(4,[1 300]));
    seq2(randperm(numel(seq2),300)) = [];
    
    % use swalign to get the alignment
    [~, al1] = seqalign(seq1,seq2);
    al1 = al1(:,1:2);
    
    p1 = al1';
    p1(1,:) = seq1(p1(1,:));
    p1(2,:) = seq2(p1(2,:));
    p1 = char(p1);
    
    % and use my C++ alignment
    al2 = swfast(seq1,seq2,[size(seq1);size(seq2)],2000);
    al2 = fillinds(al2);
    p2 = al2';
    p2(1,:) = seq1(p2(1,:));
    p2(2,:) = seq2(p2(2,:));
    p2 = char(p2);
    
    if any(sum(al1-al2 ~= 0) > 0)
        error('Indices do not agree!')
    end
end

