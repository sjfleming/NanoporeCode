function seq = RefineConsensus( consfile, bamfile, pd )

    makeall
    
    starttime = now;
    
    % now read in consensus sequence
    consseq = fastaread(consfile);
    
    % let's do baminfo/bamread for dealing with large sequences
    %bi = baminfo(bamfile,'ScanDictionary',true);
    
    ser = display_score();
    
    alignscore = blastscore('Lambda_NEB',consseq.Sequence);
    fprintf('Starting strand score: %0.1f\n',alignscore);
    ser.print(alignscore);
    
    % and get the bams
    pa = PoreAlign(bamfile,consseq.Header,consseq.Sequence);
    
    
    alparams = [];
    alparams.stripe_width = 150;
    alparams.insert_prob = 0.01;
    alparams.stay_prob = 0.02;
    alparams.extend_prob = 0.04;
    
    seq = [];

    % now loop through and refine seq a few bases at a time
    mutsize = 3000;

    
    for i=1:mutsize:numel(consseq.Sequence)
        
        % get poredata inds
        [pdinds, bminds] = pa.getInds(i:i+mutsize+500, pd);
        
        % and sort by number of matching bases, descending
        [~,indorder] = sort(-pa.Stats(bminds,2));
        pdinds = pdinds(indorder);
        bminds = bminds(indorder);
        
        % now run some refining n shit
        if numel(pdinds) > 30
            pdinds = pdinds(1:30);
        end
        
        % get the events
        events = pd.getEvents(pdinds);
        % patchy code to make deviations not matter
        for j=1:numel(events)
            foo = 0*events(j).model.sd_stdv;
            events(j).model.sd_stdv = foo + 1/sqrt(2*pi);
            events(j).model.sd_mean = foo + 0;
            events(j).stdv = 0*events(j).stdv;
        end
        % and start with the correct subsequence
        imax = min(i+mutsize+500,numel(consseq.Sequence));
        newseq = consseq.Sequence(i:imax);
            
        tic
        [~,~,events] = align_likes(newseq, events, alparams);
        fprintf('\nFull alignment took %0.1f\n',toc);
    
        for n=1:8
            seqs = viterbi_mut(newseq, events);
            [newseq,events] = MutateFast(newseq, seqs, events, alparams);
            if mod(n,2)==0
                newseq = MutateSequence(newseq, events, alparams);
            end
        end
        
        % and merge it on
        if ~isempty(seq)
            [seq,merged] = seqmerge(seq,newseq,100);
            if ~merged
                break
            end
        else
            seq = newseq;
        end
                        
        fprintf('%d bases done\n',numel(seq));
        alignscore = blastscore('Lambda_NEB',seq);
        ser.print(alignscore);
    end
    
    fprintf('%d bases total\n',numel(seq));
    alignscore = blastscore('Lambda_NEB',seq);
    ser.print(alignscore);
    
    endtime = now;
    dtime = (endtime-starttime)*24*60;
    fprintf('Assembled in %0.1f minutes\n',dtime);
    
    system(sprintf('pythonw.exe ./Python/sendmsg.py "Refined %d bases to %0.1f%% in %0.1f minutes!"',numel(seq),alignscore,dtime));
    
    % and save the output fasta file
    fid = fopen([filename '.final.fasta'],'w');
    [~,fn] = fileparts(filename);
    fprintf(fid,'>%s\n%s',fn,seq);
    fclose(fid);
    
end