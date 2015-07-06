function [score, inds] = seqalign(seq1, seq2)

    % calculate %identity
    % (with adjustments if whatsits is too big)
    [~,p,st] = swalign(seq1,seq2,'Alphabet','NT');
    
    score = 100*sum(p(2,:) == '|')/size(p,2);
    
    % and the alignment indices
    if nargout > 1
        inds = zeros([size(p,2),2]);
        inds(p(1,:) ~= '-',1) = 1:sum(p(1,:)~='-');
        inds(p(3,:) ~= '-',2) = 1:sum(p(3,:)~='-');
        inds = fillinds(inds) + repmat(st',[size(inds,1),1]) - 1;
    end
end