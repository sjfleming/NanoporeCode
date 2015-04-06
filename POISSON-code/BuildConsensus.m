function seq = BuildConsensus( filename, pd )

    makeall
    
    starttime = now;
    
    fastafile = [filename '.fasta'];
    bamfile = [fastafile '.bam'];

    % first load the fasta file, even if the data is in pd
    % (just to keep it in memory)
    seqs = fastaread(fastafile);
    
    % now let's figure out which fragments have most matched bases
    matchedbases = zeros(numel(seqs),1);
    matchedseqs = zeros(numel(seqs),1);
    
    %lambda = fastaread('.\References\Lambda_NEB.fasta');
    %lambda = lambda.Sequence;

    % now let's get some information about this here bam file
    bi = baminfo(bamfile,'NumOfReads',true,'ScanDictionary',true);
    % and scale some stuff
    matchedbases = zeros(numel(bi.ScannedDictionaryCount),1);
    for i=1:numel(bi.ScannedDictionaryCount)
        matchedbases(i) = (double(bi.ScannedDictionaryCount(i))-1)*numel(seqs(i).Sequence);
    end

	[~,inds] = sort(-matchedbases);
	seqs = seqs(inds);
    matchedbases = matchedbases(inds);
    matchedseqs = matchedseqs(inds);
    seqused = 0*matchedbases;
	
    % create PoreAlign to help with aligning etc. using first sequence
    pa = PoreAlign(bamfile,seqs(1).Header,seqs(1).Sequence);
    
    ser = display_score();
    
	% and let's build up the consensus sequence, starting with the refined
    % first sequence from PA
	seq = pa.mutateReference([],pd);
    %load('consseq.mat');
    
    alignscore = blastscore('Lambda_NEB',seq);
    fprintf('Starting strand score: %0.1f\n',alignscore);
    ser.print(100*min(numel(seq)/48502,1.0));

    
    % and an index of which reference each separate sequence came from
    seqind = 1;
    seqused(1) = 1;
    
    
    getSeq = @(h) find(strcmp({seqs.Header},h));
    
    function [newseq,newind,pa] = extendseq(side,pa)
        % find seq with enough matches, that has highest score, and also
        % aligned to left (side == 1) or right (side == 2)
        
        painds = [];
        sinds = [];
        mscores = [];
        
        for j=1:pa.NSeqs
            ind = getSeq(pa.Header{j});
            % not the same sequence, and not used?
            if strcmp(pa.Header{j},pa.RefHeader) || seqused(ind)
                continue
            end
            
            % amount overhanging to this side
            ohang = pa.Overhang(j,side);
            % and number of matches with main sequence
            nmatch = pa.Stats(j,2);
            
            if ohang < 50 || nmatch < 1000
                continue
            end
            
            painds(end+1) = j;
            sinds(end+1) = ind;
            mscores(end+1) = matchedbases(ind);
        end
        % which is der best
        
        newseq = [];
        newind = 0;
        % no best? do nothing
        if isempty(mscores)
            return
        end
        
        [ms,ind] = max(mscores);
        % index of strand in seqs
        newind = sinds(ind);
        % flag it used
        seqused(newind) = 1;
        
        isflipped = bitand(pa.Flag(painds(ind)),16);
        
        % use overhang of +400
        nbases = pa.Overhang(painds(ind),side) + 400;
        % and get the new sequence, from right or left, depending on
        % flipping setting and such, using the new porealign
        pa = PoreAlign(bamfile,seqs(newind).Header,seqs(newind).Sequence);

        % make sure it's not too long or something
        nref = numel(pa.RefSequence);
        nbases = min(nbases,nref-1);
        
        if xor(side == 1, isflipped)
            newseq = pa.mutateReference(1:nbases,pd);
        else
            newseq = pa.mutateReference(nref-nbases:nref,pd);
        end

        % now flip it if we need to before returning it
        if isflipped
            newseq = seqrcomplement(newseq);
            newind = -newind;
        end
        
        % at this point, newseq is aligned to the pa were called with
        % so the calling function needs to re-reverse it if the calling pa
        % was backwards
        
        fprintf('Added %d bases to side %d, ind = %d...\n',numel(newseq)-100,side,newind);
    end

    while 1
        
        % extend to the left side
        curind = seqind(1);

        % extend the main sequence, and simultaneously update pa
        if curind > 0
            % add a forward-aligned strand
            [newseq, newind, pa] = extendseq(1,pa);
            seqind = [newind seqind];
        else
            % add a reverse-aligned strand
            [newseq, newind, pa] = extendseq(2,pa);
            newseq = seqrcomplement(newseq);
            seqind = [-newind seqind];
        end
        
        if isempty(newseq)
            % done adding to the left
            break
        end
        
        % and merge it on
        [seq,merged] = seqmerge(seq,newseq,100);
        if ~merged
            % failed to merge, maybe just 'cause we picked the wrong
            % sequence to mutate. try again by clipping seqind
            % and resetting pa
            seqind = seqind(2:end);
            pa = PoreAlign(bamfile,seqs(abs(curind)).Header,seqs(abs(curind)).Sequence);
            fprintf('Failed to merge, trying again.\n');
            continue
        end
        
        fprintf('%d bases total\n',numel(seq));
        alignscore = blastscore('Lambda_NEB',seq);
        ser.print(100*min(numel(seq)/48502,1.0));
    end
    
    % now reset seqused, we can use the same one to extend in both
    % directions if need be
    seqused(:) = 0;
    seqused(1) = 1;
    
    % and reload porealign with first segment
    pa = PoreAlign(bamfile,seqs(1).Header,seqs(1).Sequence);
    
    while 1
        
        % extend to the right side
        curind = seqind(end);

        % extend the main sequence, and simultaneously update pa
        if curind > 0
            % add a forward-aligned strand
            [newseq, newind, pa] = extendseq(2,pa);
            seqind = [seqind newind];
        else
            % add a reverse-aligned strand
            [newseq, newind, pa] = extendseq(1,pa);
            newseq = seqrcomplement(newseq);
            seqind = [seqind -newind];
        end
        
        if isempty(newseq)
            % done adding to the right
            break
        end
        
        % and merge it on
        seq = seqmerge(seq,newseq,100);
        
        % and merge it on
        [seq,merged] = seqmerge(seq,newseq,100);
        if ~merged
            seqind = seqind(1:end-1);
            pa = PoreAlign(bamfile,seqs(abs(curind)).Header,seqs(abs(curind)).Sequence);
            fprintf('Failed to merge, trying again.\n');
            continue
        end
                
        fprintf('%d bases total\n',numel(seq));
        alignscore = blastscore('Lambda_NEB',seq);
        ser.print(100*min(numel(seq)/48502,1.0));
    end
    
    fprintf('%d bases total\n',numel(seq));
    alignscore = blastscore('Lambda_NEB',seq);
    ser.print(100*min(numel(seq)/48502,1.0));
    
    endtime = now;
    dtime = (endtime-starttime)*24*60;
    fprintf('Assembled in %0.1f minutes\n',dtime);
    
    system(sprintf('pythonw.exe ./Python/sendmsg.py "Assembled %d bases to %0.1f%% in %0.1f minutes!"',numel(seq),alignscore,dtime));
    
    % and save the output fasta file
    fid = fopen([filename '.cons.fasta'],'w');
    [~,fn] = fileparts(filename);
    fprintf(fid,'>%s\n%s',fn,seq);
    fclose(fid);
    
end