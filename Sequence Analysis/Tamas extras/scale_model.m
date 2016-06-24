function scaled_model = scale_model(model,varargin)
% SCALE_MODEL(model, options) returns a scaled version of the input model
% current levels.  'options' allows user to specify how scaling is carried
% out, and can be either:
% 'minmax' which scales model levels to data by matching the min and max
% levels, in which case this string is followed by a two element array
% containing the minimum and maximum current in the data, e.g.
%      levels = scale_model(model,'minmax',[low high]);
% 'cdf' which scales model levels to data by matching the empirical cdf of 
% each, in which case this string is followed by an array of levels, e.g.
%      levels = scale_model(model,'cdf',measured_levels);
    
    scaled_model = model;
    minmax = false;
    cdf = false;
    
    % deal with inputs
    if (numel(varargin) == 0 || numel(varargin) == 1)
        display('Not enough input arguments to ''scale_model.m''')
    elseif numel(varargin) == 2
        if strcmp(varargin{1},'minmax')
            minmax = true;
        elseif strcmp(varargin{1},'cdf')
            cdf = true;
        end
    else
        display('Too many input arguments to ''scale_model.m''')
    end
    
    % do the scaling
    if minmax
        dd = varargin{2};
        mm = [min(model) max(model)];
        scaled_model = (model-mm(1))*(dd(2)-dd(1))/(mm(2)-mm(1)) + dd(1);
    elseif cdf
        dmean = mean(varargin{2});
        dstd = std(varargin{2});
        mmean = mean(model);
        mstd = std(model);
        scaled_model = (model-mmean)*dstd/mstd + dmean;
    end

end
