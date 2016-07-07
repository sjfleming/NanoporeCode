function levels = laszlo_levels(data,p)
% laszlo_levels(data,p) implements a version of the level-finding algorithm
% documented in Laszlo et al., "Decoding long nanopore sequencing reads of
% natural DNA," Nature Biotechnology (2014), Supplementary Note 1.
% 'data' should be two columns, the first is time, the second is current
% returns a struct 'levels' containing information about each level
% Stephen Fleming, 10/7/15

    % recursive level-finding based on a probability threshold for a
    % likelihood function
    function level_search(ind1,ind2)
        
        warning('off','MATLAB:colon:nonIntegerIndex'); % turn off warning since fminbnd uses non-integer indices
        options = optimset('TolX',1);
        [possible_transition_index, min_prob] = fminbnd(@(x) log_prob(x,ind1,ind2), ind1, ind2, options);
        %probs = log_prob((ind1+1):(ind2-1),ind1,ind2);
        %[min_prob, possible_transition_index] = min(probs);
        
        if min_prob < p
            % this index is a real transition, save it and recursively look
            % at the two new levels
            level_transition_indices(a) = possible_transition_index;
            a = a+1;
            fprintf('.');
            if possible_transition_index-ind1>10
                level_search(ind1,possible_transition_index);
            end
            if ind2-possible_transition_index>10
                level_search(possible_transition_index,ind2);
            end
        end
        
    end
    
    % likelihood function giving probability that levels are the same,
    % conditioned on the index location of a possible transition
    function probability = log_prob(index,i1,i3)
        
        % with my own tweak to normalize time
        factor = 30/(size(data,1)*(data(2,1)-data(1,1)));
        probability = (index-i1) * log(nanstd(data(i1:index,2))) + ...
            (i3-index) * log(nanstd(data(index:i3,2))) - ...
            (i3-i1) * log(nanstd(data(i1:i3,2)));
        probability = probability*factor;
        
    end

% run the level search (recursive)
a = 1;
level_transition_indices = [];
level_search(1,size(data,1));

% go through each level and summarize its data, keeping it in a struct
levels = cell(numel(level_transition_indices)+1,1);
level_transition_indices = [1, round(sort(level_transition_indices)), size(data,1)];
for j = 1:(numel(level_transition_indices)-1)
    levels{j}.start_time = data(level_transition_indices(j),1);
    levels{j}.end_time = data(level_transition_indices(j+1),1);
    levels{j}.duration = data(level_transition_indices(j+1),1)-data(level_transition_indices(j),1);
    levels{j}.current_mean = nanmean(data(level_transition_indices(j):level_transition_indices(j+1),2));
    levels{j}.current_median = nanmedian(data(level_transition_indices(j):level_transition_indices(j+1),2));
    levels{j}.current_std = nanstd(data(level_transition_indices(j):level_transition_indices(j+1),2));
end

fprintf('\n')

end