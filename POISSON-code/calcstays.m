
datadir = 'C:\Minion\M13c\';

pd = PoreData(datadir);

evinds = find(pd.NumBases(:,3) > 6000 & pd.NumBases(:,3) < 9000 & pd.NumEvents(:,2) > pd.NumEvents(:,1));

alouts = zeros(numel(evinds),4);

for i=1:numel(evinds)
    % get the 2D alignment pew pew
    ald = h5read([datadir pd.Filename{evinds(i)}],'/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment');
    km = ald.kmer';
    % figure out where kmer does not advance
    km = [0; all(km(1:end-1,:) == km(2:end,:),2)];
    
    % and now total stays
    st1 = km & (ald.template>0);
    b = diff([0; st1; 0]);
    nst1 = find(b==-1) - find(b==1);
    nst1 = accumarray(nst1,1);
    try
        f = fit((1:numel(nst1))'-1,nst1,'exp1');
        alouts(i,1:2) = 100*[f.a/pd.NumEvents(evinds(i),1), exp(f.b)];
    end
    
    st2 = km & (ald.complement>0);
    b = diff([0; st2; 0]);
    nst2 = find(b==-1) - find(b==1);
    nst2 = accumarray(nst2,1);
    try
        f = fit((1:numel(nst2))'-1,nst2,'exp1');
        alouts(i,3:4) = 100*[f.a/pd.NumEvents(evinds(i),2), exp(f.b)];
    end
    
    i
end

evouts = 0*alouts;
evouts(:,1) = alouts(:,1)./(1-0.01*alouts(:,2));
evouts(:,2) = alouts(:,3)./(1-0.01*alouts(:,4));
evouts(:,3:4) = 100*pd.NumStays(evinds,:)./pd.NumEvents(evinds,:);

alouts = alouts(all(abs(alouts) > 1,2),:);
mean(alouts)
std(alouts)
