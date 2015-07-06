function score = blastscore(refname, seq)
% use blastn to blast-compare two sequences

    fn1 = ['.\References\' refname '.fasta'];
    fn2 = '.\seq.fasta';
    
    % save sequence to a temporary file
    fid = fopen(fn2,'w');
    fprintf(fid,'>Sequence\n');
    fprintf(fid,seq);
    fclose(fid);
    
    % -q mismatch, -r match, -G gap start, -E gap extend
    
    % execute blast program
    runstr = sprintf('.\\Tools\\bl2seq.exe -r 1 -q -2 -G 5 -E 2 -i %s -j .\\seq.fasta -p blastn',fn1);
    %runstr = sprintf('.\\Tools\\blastn.exe -reward 1 -penalty -2 -gapopen 5 -gapextend 2 -evalue 0.1 -query %s -subject .\\seq.fasta',fn1);
    [~, cmdout] = system(runstr);
    
    % find all identity scores from output string
    [~,stok] = regexp(cmdout, 'Identities = (\d+)/(\d+)', 'match', 'tokens');
    % and the expectation values
    [~,etok] = regexp(cmdout, 'Expect = (\d+\.\d+)', 'match', 'tokens');
    
    % now remove temporary file
    delete(fn2);
    
    if isempty(stok)
        fprintf(2,'No blast matches found!\n');
        score = 0;
        return
    end
    
    % now go through and add up everything with low e-values
    nmatch = 0;
    ntot = 0;
    for i=1:numel(stok)
        escore = str2double(etok{i}{1});
        if escore > 0.1
            break;
        end
        nmatch = nmatch + str2double(stok{i}{1});
        ntot = ntot + str2double(stok{i}{2});
    end
    
    % and return the best match score
    score = 100*nmatch/ntot;
    
    % and display a cute little thing
    fprintf(2,'Blast score: %0.1f with %d bases ------------\n',score,ntot);

end