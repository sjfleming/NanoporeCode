function [ filtdata ] = filt_decimate( data, r )
    %FILT_DECIMATE Filters data using downsampling decimation
    %   Factor r is number of original points per downsampled point.
    %   Data is passed with colums of [time, sig1, sig2, ...], and
    %       returns the same array but with the data filtered.

    if size(data,2)>10
        data = data';
    end
    
    for i = 1:size(data,2) % go through each column of data separately
        d = data(:,i);
        filtdata(:,i) = accumarray(1+floor((1:numel(d))/r)',d',[],@mean);
%         
%         % Matlab says not to decimate with r > 13, instead do it multiple times
%         reps = floor(log10(r));
%         left = ceil(r/(10^reps));
%         if r < 13
%             f = decimate(data(:,i),r,'fir');
%         elseif reps>1
%             f = decimate(data(:,i),10,'fir');
%             for j=2:reps
%                 if numel(f)>100
%                     f = decimate(f,10,'fir');
%                 end
%             end
%             if left > 1
%                 f = accumarray(1+floor((1:numel(f))/left)',f',[],@median);
%             end
%         else
%             f = decimate(data(:,i),10,'fir');
%             if left > 1
%                 f = decimate(f,left,'fir');
%             end
%         end
%         
%         filtdata(:,i) = f';
%         
     end
    
end
