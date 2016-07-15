function levels = find_levels_ks(mol, filter, sample, p, tr)
% find_levels_ks(data,p) implements a recursive level-finding algorithm
% based on Kolmogorov-Smirnov testing.
% p-value should be 0.07, determined empirically.

% mol is a molecule object.  it should have the start and end times filled in.
% filter is the filter frequency
% sample is the sampling frequency
% these are used for downsampling the original signal

% returns a struct 'levels' containing information about each level
% Stephen Fleming, 7/1/16

% NOTE: doesn't work starting from one big level, because randomly sampling
% from a huge bunch of levels will never lead to convincingly different
% CDFs, they will just be all blurry and never cross the threshhold...

    % recursive level-finding based on a probability threshold for a
    % likelihood function
    function level_search(ind1,ind2)
        
        warning('off','MATLAB:colon:nonIntegerIndex'); % turn off warning since fminbnd uses non-integer indices
        options = optimset('TolX',1);
        
%         figure(8)
%         hold on
%         xx = round(ind1)+20:round(ind2)-20;
%         yy = arrayfun(@(x) inv_ks_stat(x,ind1,ind2), xx);
%         plot(xx,yy)
        
        % find index which minimizes the probability that levels are the same
        [possible_transition_index, min_inv_ks_stat] = fminbnd(@(x) inv_ks_stat(x,ind1,ind2), ind1, ind2, options);
        
        if min_inv_ks_stat < p % empirically determined cutoff
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
    function kinv = inv_ks_stat(index,i1,i3)
        
        d1 = randsample(data(i1:index,2),floor(min(2000,index-i1+1)));
        d2 = randsample(data(index:i3,2),floor(min(2000,i3-index+1)));
        [~,~,k] = kstest2(d1,d2);
        kinv = 1-k;
        
    end

% do chunks if necessary, so we can sample the data at 5kHz
f = util.getMoleculeFilesAndTimes(mol);
totaltime = sum(f(:,3)-f(:,2));
numpts = sample*totaltime; % total # pts needed, sampled at 5kHz
chunks = ceil(numpts / 2^22); % number of max 2^22 point data chunks necessary
% times = linspace(0, totaltime, chunks+1)';
% times(1:end-1,2) = times(2:end,1);
% times = times(1:end-1,:); % row matrix of chunks of time to grab
% % times = [55.8, 58];
times = tr;

% for each chunk of data
for i = 1:size(times,1)
    
    % grab the data from the molecule
    trange = times(i,:);
    data = util.doLoadMoleculeData(mol, diff(trange)*sample, 'minmax', filter, trange); % 2 kHz filter, minmax
    
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

end