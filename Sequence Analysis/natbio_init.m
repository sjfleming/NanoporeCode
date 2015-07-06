function [events,M13] = natbio_init(seed, count)

    make all
    
    M13 = fastaread('./References/M13mp18_EcoRI.fasta');
    M13 = M13.Sequence;
    
    if nargin < 1
        events = [];
        return
    end

    if ~isunix()
        pd = PoreData('C:\Minion\M13c');
    else
        pd = PoreData('~/Minion/M13c');
    end

    
    % now we want to find all "useable" strands
    % these are defined as having a reasonable number of levels
    
    evinds = find(pd.NumBases(:,3) > 6000 & pd.NumBases(:,3) < 8000 & pd.NumEvents(:,2) > pd.NumEvents(:,1));
    % and randomize them, but not too randomly
    rng(seed);
    evinds = evinds(randperm(numel(evinds)));
    
    % get a subset
    curinds = evinds(1:count);
    
    % and actually get the events
    events = pd.getEvents(curinds);

end