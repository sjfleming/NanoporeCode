function make(arg)

    basedir = './cpp/';
    ext = ['.' mexext];

    function time = filetime(file)
        finfo = dir(file);
        if isempty(finfo)
            time = 0;
        else
            time = finfo.datenum;
        end
    end

    function time = recents(files)
        time = 0;
        for j=1:numel(files)
            time = max(time,filetime([basedir files{j}]));
        end
    end

    function domake(srcs, deps)
        time = recents(strsplit(srcs));
        time = max(time,recents(strsplit(deps)));
        outfile = strsplit(srcs);
        outfile = strsplit(outfile{1},'.');
        outfile = [outfile{1} ext];

        if time > filetime(outfile)
            srcs = strsplit(srcs);
            for j=1:numel(srcs)
                srcs{j} = [basedir srcs{j}];
            end
            % make the thing
            mex('CXXFLAGS=$CXXFLAGS -std=c++0x',srcs{:});
        end
    end

    switch (arg)

        case 'all'
            % Compile all changed files and whatnot

            %mexfiles = {'viterbi1d','viterbi2d','viterbikd','align_like','align_fast','align_local','align_like_omg'};
            alignh = 'Alignment.h EventData.h AlignUtil.h Sequence.h MutateAlign.h';
            domake('align_likes.cpp Alignment.cpp hat.cpp',alignh)
            domake('MutateSequence.cpp MutateAlign.cpp Alignment.cpp hat.cpp',alignh)
            domake('ScoreMutations.cpp MutateAlign.cpp Alignment.cpp hat.cpp',alignh)
            domake('viterbikd.cpp',alignh)
            domake('cigarstats.cpp','')
            domake('cigaroverlap.cpp','')
            domake('swfast.cpp','')

        case 'clean'
            % Remove all changed files
            delete(['./*.' mexext])
    end
end
