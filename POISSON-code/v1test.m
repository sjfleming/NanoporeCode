function p = v1test() 

    pd = PoreData('C:\Minion\Lambda-burnin\');

    evind = 22;

    evt = pd.getEvent(evind,'t');
    evt.ref_align = (1:numel(evt.ref_align))';
    evc = pd.getEvent(evind,'c');
    strt = pd.getSequence(evind,'t');
    strc = seqrcomplement(pd.getSequence(evind,'c'));
    str = pd.getSequence(evind,'2d');

    mex viterbikd.cpp

    tic
    params = [];
    params.mutations = 0;
    params.skip_prob = 0.10;
    params.stay_prob = 0.01;
    dpath = viterbikd(evt, params);
    toc

    seq = statestoseq(dpath(:,1));
    [s,p] = nwalign(seq,strt,'Alphabet','NT')
    
end