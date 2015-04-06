function [ filtdata ] = filt_decimate( data, r )
    %FILT_DECIMATE Filters data using downsampling decimation
    %   Factor r is number of original points per downsampled point.
    %   Data is passed with colums of [time, sig1, sig2, ...], and
    %       returns the same array but with the data filtered.

    %filtdata = accumarray(1+floor((1:numel(data))/r)',data',[],@median)';
    
    % Matlab says not to decimate with r > 13, instead do it multiple times
    reps = round(r/10);
    left = round(r/(10*reps));
    if r < 13
        filtdata = decimate(data,r);
    else
        filtdata = decimate(data,10);
        for i=2:reps
            filtdata = decimate(filtdata,10);
        end
        if left > 1
            filtdata = decimate(filtdata,left);
        end
    end
    
end
