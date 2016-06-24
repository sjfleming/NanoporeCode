function al2test( )

    lambda = fastaread('.\References\Lambda_NEB.fasta');
    lambda = lambda.Sequence;

    pd = PoreData('C:\Minion\Lambda-burnin\');
    
    evind = 19;
    
    evt = pd.getEvent(evind,'t');
    evc = pd.getEvent(evind,'c');
    % get the sequences
    strt = pd.getSequence(evind,'t');
    % reversed
    strc = seqrcomplement(pd.getSequence(evind,'c'));
    str = pd.getSequence(evind,'2d');
    
    % find out where they go
    [lstr,~,match] = seqintersect(lambda, str);
    
    mex align_like.cpp
    mex viterbikd.cpp
    
    events = [evt evc];
    
    events(1).ref_align = (1:numel(events(1).ref_align))';
    events(2).ref_align = 0*events(2).ref_align;

    
    function [s, newev] = testseq(states)
        [s1, newev(1)] = align_like(events(1),states);
        [s2, newev(2)] = align_like(events(2),states);
        fprintf('Aligned (%d, %d)\n',sum(newev(1).ref_align>0),sum(newev(2).ref_align>0));
        s = s1 + s2;
    end
    
    % now try a few things, see how their scores match up
    %testseq(strt)
    %testseq(strc)
    %testseq(str)
    %testseq(lstr)
            
    for n=1:50
        dpaths = viterbikd(events);
        events(1).ref_align = 0*events(1).ref_align;
        events(2).ref_align = 0*events(2).ref_align;

        scores = zeros(size(dpaths,2),1);
        newevents = {};
        for i=1:size(dpaths,2)
            [scores(i), newevents{i}] = testseq(dpaths(:,i));
        end
        [m,ind] = max(scores(2:end));
        ind = ind + 1;
        fprintf('Best: %0.1f | Max: %0.1f at %d\n',scores(1),m,ind);
        seq = statestoseq(dpaths(:,ind));
        states = seqtostates(seq);
        
        fprintf('Best: %0.1f%% | Using: %0.1f%%\n',seqalign(statestoseq(dpaths(:,1)),lstr),...
                seqalign(statestoseq(states),lstr));
            
        seq1 = mutateseqs(statestoseq(dpaths(:,1)),[newevents{1}(1) newevents{ind}(1)],dpaths(:,[1 ind]));

        [~,events(1)] = align_like(events(1),states);
        [~,events(2)] = align_like(events(2),states);
        
        % and fill ref_align back in, don't want to throw away events
        [m,ind] = min(events(1).ref_align);
        while (ind > 0 && m > 0)
            events(1).ref_align(ind) = m;
            ind = ind - 1;
            m = m - 1;
        end
        [m,ind] = max(events(1).ref_align);
        while ind <= numel(events(1).mean)
            events(1).ref_align(ind) = m;
            ind = ind + 1;
            m = m + 1;
        end
    end
end