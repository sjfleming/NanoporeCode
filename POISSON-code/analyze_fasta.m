function analyze_fasta(pathname)
%ANALYZE_FASTA outputs fasta files from fast5 files
    % trim path separator from end of path
    if pathname(end) == '/' || pathname(end) == '\'
        pathname = pathname(1:end-1);
    end
    
    % get file metadata
    pd = PoreData(pathname);
    
    % locations for stuff
    locs = {'template','complement','2D'};
    
    fastapath = cell(1,3); % filename
    fastafid = cell(1,3);  % file id
    
    % make directories and paths and files and stuff
    for i=1:3
        fastapath{i} = [pathname '_' locs{i} '.fasta'];
        fastafid{i} = fopen(fastapath{i},'w');
    end
        
    tic;
    for i=1:pd.NumFiles
        
        for j=1:3
            fwrite(fastafid{j}, pd.getFasta(i,locs{j}));
        end
        
        if mod(i,100) == 0
            t = toc;
            trem = (pd.NumFiles-i)*t/i;
            fprintf(1,'%0.2f minutes remaining...\n',trem/60);
        end
    end
    
    for i=1:3
        fclose(fastafid{i});
    end

end

