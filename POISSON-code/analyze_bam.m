function idents = analyze_bam(pathname)
%ANALYZE_BAM calculate some statistics on BAM files

    tic
    
    % check for .bam-suffixed files
    locs = {'template','complement','2D'};
    
    % trim path separator from end of path
    if pathname(end) == '/' || pathname(end) == '\'
        pathname = pathname(1:end-1);
    end
    
    [~,filebase,~] = fileparts(pathname);
    
    idents = cell(numel(locs),1);
    
    for j=1:numel(locs)
        
        filename = [pathname '_' locs{j} '.fasta.bam'];

        % first, load bam stuff into memory
        try
            h=baminfo(filename,'ScanDictionary',true);
        catch
            fprintf(2,'Not found: %s\n',filename);
            continue
        end
        bm = BioMap(filename, 'SelectReference', h.ScannedDictionary{1});

        % get coverage info
        cover = bm.getBaseCoverage(1,h.SequenceDictionary.SequenceLength);
        fprintf('Aligned reads (%s): %d\n',locs{j},bm.NSeqs);
        fprintf('Mean coverage (%s): %0.2f\n',locs{j},mean(cover));

        % get the reference sequence
        try
            refseq = fastaread([cd() '\\References\\' h.ScannedDictionary{1} '.fasta']);
        catch
            error('Reference file not found!')
        end
        refseq = refseq.Sequence;

        % now loop through and calculate some stuff
        ids = [];
        
        for i=1:bm.NSeqs

            % calculate statistics from cigar string
            st = cigarstats(bm.Signature{i},refseq,bm.Sequence{i},bm.Start(i));
            
            % and how many matches we have/how many in ref (% identity)
            % (this is basically accuracy)
            id.accuracy = 100*st(2)/st(1);
            
            % the name
            id.name = bm.Header{i};
            
            % start and end indices
            id.start = bm.Start(i);
            id.end = bm.Start(i) + sum(st([2 3 5]));
        
            if numel(ids)==0
                ids = id;
            else
                ids(end+1) = id;
            end
        end
        
        idents{j} = ids;
    end
    
    toc

    % make nice plot
    xspc = 3.0;
    xs = xspc/2:xspc:100-xspc/2;
    
    nplots = sum(1-cellfun(@isempty,idents));
    
    if nplots==0
        return
    end
    
    delts = [0 0 0];
    if nplots == 2
        delts = [-xspc/4 0 xspc/4];
    elseif nplots == 3
        delts = [-xspc/3 0 xspc/3];
    end
    
    colors = [194 224 250; 210 226 139; 255 154 200]/255.0;
    legs = {};
    
    fig = figure;
    set(fig, 'Position', [600 290 641 586]);
    
    for j=1:3
        if isempty(idents{j})
            continue
        end
        
        counts = hist([idents{j}.accuracy], xs);
        bar(xs+delts(j),counts,1.0/nplots,'FaceColor',colors(j,:))
        hold on
        legs{end+1} = locs{j};
    end
    hold off
    ylabel('# of aligned strands','FontSize',13);
    xlabel('% identity','FontSize',13);
    title([filebase ' alignment results'],'FontSize',13);
    legend(legs)
    xlim([50 100])
    options.Format = 'png';
    
    hgexport(gcf,[pathname '.png'],options);
    
end

