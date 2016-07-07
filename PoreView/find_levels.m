function levels = find_levels(mol, filter, sample, p)
% find_levels(data,p) implements a version of the level-finding algorithm
% documented in Laszlo et al., "Decoding long nanopore sequencing reads of
% natural DNA," Nature Biotechnology (2014), Supplementary Note 1.
% BUT, it is altered by me to do a better job.

% mol is a molecule object.  it should have the start and end filled in.

% 'data' should be two columns, the first time, the second current
% returns a struct 'levels' containing information about each level
% Stephen Fleming, 6/28/16

    % recursive level-finding based on a probability threshold for a
    % likelihood function
    function level_search(ind1,ind2)
        
        warning('off','MATLAB:colon:nonIntegerIndex'); % turn off warning since fminbnd uses non-integer indices
        options = optimset('TolX',1);
        
%         figure(8)
%         hold on
%         xx = round(ind1)+20:round(ind2)-20;
%         yy = arrayfun(@(x) log_prob(x,ind1,ind2), xx);
%         plot(xx,yy)
        
        % find index which minimizes the probability that levels are the same
        [possible_transition_index, min_prob] = fminbnd(@(x) log_prob(x,ind1,ind2), ind1, ind2, options);
        
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
        
        sig = nanstd(data(i1:i3,3));
        sig1 = nanstd(data(i1:index,3));
        sig2 = nanstd(data(index:i3,3));
        mu = nanmean(data(i1:i3,2));
        mu1 = nanmean(data(i1:index,2));
        mu2 = nanmean(data(index:i3,2));
        
        % probability that all of these points are from two Gaussians
        % versus one Gaussian
        laszlo_updated = ((index-i1) * sig1 + ...
            (i3-index) * sig2 - ...
            (i3-i1) * sig) / (2*sig^2) - log(sig)/2 - log(2*pi);
        
        % Welch's t-test
        t = (mu1 - mu2) / sqrt(sig1^2/(index-i1) + sig2^2/(i3-index));
        welch = log(tcdf(abs(t),3,'upper')); % hypothesis testing with t-distribution, 3 degrees of freedom
        
        % probability that the resulting distributions actually have the
        % same mean
%         fleming = -(mu1 - mu2)^2  * ...
%             (sig1^2 + sig2^2) / ...
%             (2 * sig1^2 * sig2^2);
        fleming = (-(mu1-mu)^2 - (mu2-mu)^2) / (2 * sig^2) * (index-i1)*(i3-index)/(i3-i1)^2;
        
        probability = min(0,laszlo_updated) + min(0,welch) + min(0,fleming);
        
    end

% do chunks if necessary, so we can sample the data at 5kHz
f = util.getMoleculeFilesAndTimes(mol);
totaltime = sum(f(:,3)-f(:,2));
numpts = sample*totaltime; % total # pts needed, sampled at 5kHz
chunks = ceil(numpts / 2^22); % number of max 2^22 point data chunks necessary
times = linspace(0, totaltime, chunks+1)';
times(1:end-1,2) = times(2:end,1);
times = times(1:end-1,:); % row matrix of chunks of time to grab
% times = [55.8, 58];

% for each chunk of data
for i = 1:size(times,1)
    
    % grab the data from the molecule
    trange = times(i,:);
    data = util.doLoadMoleculeData(mol, diff(trange)*sample, 'pointwise', 500, trange); % 500 Hz filter
    data_filt = util.doLoadMoleculeData(mol, diff(trange)*sample, 'minmax', filter, trange); % 2 kHz filter, minmax
    data(:,3) = data_filt(1:size(data,1),2);
    clear data_filt;
    
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
        levels{j}.current_mean = nanmean(data(level_transition_indices(j):level_transition_indices(j+1),3));
        levels{j}.current_median = nanmedian(data(level_transition_indices(j):level_transition_indices(j+1),3));
        levels{j}.current_std = nanstd(data(level_transition_indices(j):level_transition_indices(j+1),3));
    end
    
    fprintf('\n')
    
end

end