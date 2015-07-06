function align_fasta(fastabase, refname)

    locs = {'template','complement','2d'};
    
    curdir = cd;
    basedir = [curdir '\'];

    % go to cygwin dir for tools dependencies
    cd('c:\cygwin\bin');
    
    cdsp = strsplit(curdir,{'/','\'});
    cdsp{1} = '/cygdrive/c';
    cdsp{end+1} = '';
    % first, turn fastapath into cygwin path
    cygbase = strjoin(cdsp,'/');


    % trim path separator from end of path
    if fastabase(end) == '/' || fastabase(end) == '\'
        fastabase = fastabase(1:end-1);
    end

    fastasp = strsplit(fastabase,{'/','\'});
    fastasp{1} = '/cygdrive/c';

    % first, turn fastapath into cygwin path
    fastabase = strjoin(fastasp,'/');
    
    for i=1:3
        
        fastapath = [fastabase '_' locs{i} '.fasta'];

        reffile = [cygbase 'References/' refname '.fasta'];
        indfile = [cygbase 'References/' refname '-index'];

        samtools = [basedir 'Tools\\samtools.exe'];

        % build indexed version of reference
        system([samtools ' faidx ' reffile]);
        % and a database version
        system([basedir 'Tools\\lastdb.exe ' indfile ' ' reffile]);
        fprintf('Reference indexed, starting alignment...\n');
        % and then do the alignment
        runstr = [basedir 'Tools\\lastal.exe -s 2 -T 0 -Q 0 -a 1 ' indfile ' ' fastapath ' | '...
            'python.exe ' basedir 'Tools\\maf-convert.py sam | ' samtools ' view -t ' reffile '.fai -S -b - | '...
            samtools ' sort - ' fastapath];
        tic
        system(runstr);
        toc
        % and finally index the outputted file
        system([samtools ' index ' fastapath '.bam']);
        fprintf('Alignment done and indexed.\n');
        
    end

    cd(curdir);
    
end

