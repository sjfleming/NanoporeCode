function seq=v2test()

    lambda = fastaread('.\References\Lambda_NEB.fasta');
    lambda = lambda.Sequence;

    pd = PoreData('C:\Minion\Lambda-burnin\');
    
    evind = 22;
    
    evt = pd.getEvent(evind,'t');
    evc = pd.getEvent(evind,'c');
    strt = pd.getSequence(evind,'t');
    strc = seqrcomplement(pd.getSequence(evind,'c'));
    str = pd.getSequence(evind,'2d');
    
    al = pd.getAlignment(evind);
    
    mex viterbi2d.cpp

    tic
    [scores,dpath]=viterbi2d(evt,evc,mod1,mod2);
    toc
    
    % now create function to respond to clickies with dpaths
    
    scores = exp(scores);
    pct = prctile(rshape(scores),[5,95]);
    
    fig = figure(1);
    h = imagesc(scores,'HitTest','off');
    ax = get(h,'Parent');
    
    foo = [];
    foo.scores = scores;
    foo.dpath = dpath;
    set(ax,'UserData',foo);
    
    set(ax,'CLim',pct);
    
    colormap gray
    hold on

    plot(dpath(:,2)+1,dpath(:,1)+1,'r','HitTest','off','Tag','dpath')
    plot(al(:,2),al(:,1),'b')
    plot([0 pd.NumEvents(evind,2)],[0 pd.NumEvents(evind,1)],'g')
    
    seq = statestoseq(dpath(:,3))
    [s,p]= nwalign(str(1:numel(seq)),seq,'Alphabet','NT');
    s
    p
    [s,p]= swalign(seqrcomplement(lambda),seq,'Alphabet','NT');
    s
    p
    [s,p]= swalign(seqrcomplement(lambda),str(1:numel(seq)),'Alphabet','NT');
    s
    p
end