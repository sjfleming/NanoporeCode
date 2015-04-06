function seq = statestoseq(states)
% Convert a sequence of states to the original sequence

    if size(states,2) == 1
        seq = statestoseq1(states);
    else
        % loop through
        seq = {};
        for i=1:size(states,2)
            seq{end+1} = statestoseq1(states(:,i));
        end
    end
end

% internal function that just returns one sequence
function seq = statestoseq1(states)
nstates = 1024;
    
    % states come in 1-based, so subtract 1
    states = states - 1;

    % states are zero-based (0...1023) for this part
    function st=next_state(state, ind, nsteps)
        st = bitand(bitshift(state, 2*nsteps),nstates-1)+ind;
    end
    % return base at ind, ind = 0...4, 0 is leftmost, 4 is rightmost
    function base=get_base(state, ind)
        base = bitand(3,bitshift(state,2*(ind-4)));
    end

    % first, get them out numerically, then convert to letters
    % start by getting the first part of the first kmer
    seq = [];
    
    curstate = states(1);
    seq(1) = get_base(states(1),0);

    
    nmismatch = 0;
    nskip = 0;
    nstay = 0;
    
    % now loop through the rest of the states
    for i=2:numel(states)
        if curstate == states(i)
            % we stayed where we were, no transition
            % curstate stays, states stays, nothing added to sequence
            nstay = nstay + 1;
            continue
        end
        % otherwise, we have a new state, figure out what the possible
        % transition could have been
        for nskips=1:4
            for ind=0:4^nskips
                if next_state(curstate,ind,nskips)==states(i)
                    % we found it! first, insert missing letters
                    % already added curstate's zeroth base
                    for j=1:nskips
                        seq(end+1) = get_base(curstate,j);
                    end
                    % now update curstate
                    curstate = states(i);
                    % and write relevant vars, '1' skip is 0
                    % i really apologize for this line of code
                    % the naming and indexing both suck
                    nskip = nskip + nskips - 1;
                    break
                end
            end
            if curstate==states(i)
                % we found it in inner loop
                break
            end
        end
        % didn't find it whaaa?
        if curstate ~= states(i)
            % likely we have a mismatch
            % fuck it, we'll do it live
            curstate = states(i);
            % take the middle letter of it
            seq(end+1) = get_base(curstate,0);
            nmismatch = nmismatch+1;
        end
    end
    
    if (nmismatch > 0)
        fprintf(2,'states->seq with %d mismatches!\n',nmismatch);
    end
    %fprintf(1,'states->seq with %d skips and %d stays\n',nskip,nstay);
    
    % now pop on the last bases we didn't do yet
    for i=1:4
        seq(end+1) = get_base(curstate,i);
    end

    % convert it to char/sequence array
    seq(seq==0) = 'A';
    seq(seq==1) = 'C';
    seq(seq==2) = 'G';
    seq(seq==3) = 'T';
    seq = char(seq);
end