classdef PoreAlign
    %POREALIGN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        RefSequence     % sequence of reference
        RefHeader       % header name of reference
        
        BamFile         % filename
        
        Overhang        % how much each strand goes over end of ref
        Stats           % cigar stats, saved
        
        NSeqs
        Signature
        Sequence
        Start
        End
        Flag
        Header
        
    end
        
    properties (Access = private)
        
        bm              % BioMap object
        
    end
    
    methods
        
        function obj = PoreAlign(bamfile, refname, refseq)
            % Build BioMap superclass
            
            if nargin < 2
                h=baminfo(bamfile,'ScanDictionary',true);
                refname = h.ScannedDictionary{1};
            end
            if nargin < 3
                refseq = [];
            end
            
            obj.bm = BioMap(bamfile, 'SelectReference', refname, 'InMemory', true);
            
            % set copies of bm
            obj.NSeqs = obj.bm.NSeqs;
            obj.Signature = obj.bm.Signature;
            obj.Sequence = obj.bm.Sequence;
            obj.Start = obj.bm.Start;
            obj.Flag = obj.bm.Flag;
            obj.Header = obj.bm.Header;

            obj.RefSequence = refseq;
            obj.RefHeader = refname;
            obj.BamFile = bamfile;
            
            % calculate overhangs, ends, etc.
            obj.Overhang = zeros(obj.NSeqs,2);
            obj.Stats = zeros(obj.NSeqs,7);
            obj.End = zeros(obj.NSeqs,1);
            for i=1:obj.NSeqs
                st = cigarstats(obj.Signature{i},obj.RefSequence,obj.Sequence{i},obj.Start(i));
                obj.Stats(i,:) = st;
                % the end index is simply the number in the reference
                % sequence that were traversed during the cigar
                % aka start+match+mismatch+deletions
                obj.End(i) = obj.Start(i) + sum(st([2 3 5]));
                % then the overhangs are just start minus hard clips, more
                % or less (clipped to 0)
                obj.Overhang(i,:) = max(0,[st(6)-obj.Start(i),st(7)-numel(obj.RefSequence)+obj.End(i)]);
            end
            
        end
        
        function disp(obj)
            % Overload of internal display method
            fprintf('PoreAlign\n\n\t\t%s\n\t\t%s\n\t\t%d aligned reads\n',obj.BamFile,obj.RefHeader,obj.bm.NSeqs);
        end
        
        function [pdinds, bminds] = getInds(obj, seqrange, pd)
            % Find indices that align with a range of reference, return
            % with correct signs
            
            seqrange = [min(seqrange),max(seqrange)];
            
            % get the indices
            bminds = obj.bm.getIndex(seqrange(1),seqrange(2),'overlap',diff(seqrange)/2);
            
            % first, remove the one that is equal to itself
            bminds = bminds(~strcmp(obj.Header(bminds),obj.RefHeader));
            
            % get the corresponding poredata indices
            pdinds = pd.getIndex(obj.Header(bminds));
            
            % and flip the negatives
            pdinds(obj.Flag(bminds) == 16) = -pdinds(obj.Flag(bminds) == 16);
            
        end
        
        function seq = mutateReference(obj, seqrange, pd)
            % Improves reference fragment via mutation
            
            % condition inputs
            if isempty(seqrange)
                seqrange = [1 numel(obj.RefSequence)];
            end

            seqrange = [min(seqrange),max(seqrange)];
            
            seqrange(1) = max(seqrange(1),1);
            seqrange(2) = min(seqrange(2),numel(obj.RefSequence));

            lambda = fastaread('.\References\Lambda_NEB.fasta');
            lambda = lambda.Sequence;

            mutsize = 3000;
            
            % are we biting off more than we can chew?
            if diff(seqrange) > mutsize + 500
                % do a shorter piece first, and extend it a bit to get overlap
                seq1 = obj.mutateReference([seqrange(1),seqrange(1)+mutsize+400],pd);
                % then the longer piece, which'll probably subdivide
                seq2 = obj.mutateReference([seqrange(1)+mutsize,seqrange(2)],pd);
                % and then merge them
                seq = seqmerge(seq2,seq1);
                % and return here
                return
            end
            
            fprintf('Mutating reference from %d to %d\n',seqrange(1),seqrange(2));
            
            % get PoreData indices of matching events, incl. negatives
            [pdinds, bminds] = obj.getInds(seqrange, pd);
            
            
            % use our cigaroverlap function to get matching subsequence
            seqs = {};
            seqns = [];
            for i=1:numel(bminds)
                seqs{i} = cigaroverlap(obj.Signature{bminds(i)},obj.Sequence{bminds(i)},obj.Start(bminds(i)),seqrange(1),seqrange(2));
                seqns(i) = sum(seqs{i} ~= '-');
            end
            
            % keep 12 seqs with length closest to target/ref
            %[~,inds] = sort(abs(diff(seqrange)-seqns));
            
            fprintf('Found %d aligned strands\n',numel(pdinds));
            % keep the ones with the best oxford scores
            if numel(pdinds) > 20
                scores = pd.Score(abs(pdinds),3);
                [~,inds] = sort(-scores);
                pdinds = pdinds(inds);
                pdinds = pdinds(1:20);
            end
                        
            % get the events
            events = pd.getEvents(pdinds);
            % and start with the correct subsequence
            seq = obj.RefSequence(seqrange(1):seqrange(2));
            
            fprintf('Initial score: %0.1f\n',seqalign(lambda,seq));
            
            fprintf('Doing full alignment with %d...\n',numel(events));
            
            alparams = [];
            alparams.stripe_width = 150;
            alparams.insert_prob = 0.01;
            alparams.stay_prob = 0.02;
            alparams.extend_prob = 0.04;
            alparams.do_fast = true;

            tic
            events = seedaligns(events, seq, alparams);
            fprintf('\nSeeding took %0.1f\n',toc);
    
            % stuff.
            for n=1:5
                % take some random events
                [ev, evinds] = randsubset(events,20);
                ss = randsubset(seqs,10);
                [seq,ev] = MutateFast(seq, ss, ev, alparams);
                %MutateTest(seq,ev,alparams);
                events(evinds) = ev;
            end
            % and display the score for this fragment
            % shouldn't really break anything I hope
            fprintf('Final score: %0.1f\n',seqalign(lambda,seq));
            
        end

    end
    
end