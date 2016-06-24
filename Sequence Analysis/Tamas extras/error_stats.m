function [homo,ins,del,mut] = error_stats(seq1, seq2)

    [~,p,st] = swalign(seq1,seq2,'Alphabet','NT');
    
    % first, find locations of homopolymer regions in seq1
    bases = 'ACGT';
    hpinds = [];
    for i=1:4
        hpinds = [hpinds strfind(seq1,repmat(bases(i),[1,5]))];
    end
    hpinds = sort(hpinds);
    hpstr = seq1;
    hpstr(:) = '0';
    for i=1:numel(hpinds)
        hpstr(hpinds(i):hpinds(i)+4)='1';
    end
    % consider near HP as an HP?
    %nearhp = imdilate(hpstr == '1',[1 1 1]);
    %hpstr(nearhp) = '1';
    % now find which errors correspond to homopolymer errors
    % by substituting into p
    pnew = p;
    pnew(1,pnew(1,:) ~= '-') = hpstr(st(1):st(1)+sum(pnew(1,:) ~= '-')-1);
    % and count
    errinds = find(p(2,:) ~= '|');
    % near homopolymer region
    % (deletions only)
    homo = accumarray(nt2int(p(1,p(3,:)=='-' & pnew(1,:)=='1'))',1,[4 1]);
    % insertion
    ins = accumarray(nt2int(p(3,p(1,:)=='-'))',1,[4 1]);
    % deletions
    del = accumarray(nt2int(p(1,p(3,:)=='-'))',1,[4 1]);
    % non-homopolymer-related only
    del = del - homo;
    % mutations
    p([1,3],all(p([1,3],errinds)~='-'));
    mut = accumarray(nt2int(p([1,3],errinds(all(p([1,3],errinds)~='-'))))',1,[4 4]);
    
    fprintf('%d errors total\n',numel(errinds));
    fprintf('%d homopolymer\n',sum(homo));
    fprintf('%d insertions\n',sum(ins));
    fprintf('%d deletions\n',sum(del));
    fprintf('%d mismatches\n',sum(sum(mut)));
end