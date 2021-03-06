function [ filtdata ] = filt_lpb( data, n, fp )
    %FILT_LPB Filters data using a low-pass Bessel filter
    %   Frequency fp is in Hz, n is poles.
    %   Data is passed with colums of [time, sig1, sig2, ...], and
    %       returns the same array but with the data filtered.


    % extract si. easier than passing each time
    si = data(2,1)-data(1,1);
    fs = 1/si;

    % create the (analog) filter coefficients
    [z,p,k] = besself(n, 2*pi*fp);
    % convert to digital filter zero-pole-gain
    [zd,pd,kd]=bilinear(z,p,k,fs,fp);
    % and create second-order filter cascade for calculating
    %sos=real(zp2sos(zd,pd,kd));
    sos=zp2sos(zd,pd,kd);

    % copy
    filtdata = data;
    % and replace data portion for all dimensions
    for i=2:size(filtdata,2)
        filtdata(:,i) = real(filtfilt(sos,1,data(:,i)));
    end


%     % extract si. easier than passing each time
%     si = data(2,1)-data(1,1);
% 
%     % convert from absolute frequency to 'normalized frequency'
%     % which is 0 at 0 and 1 at Nyquist freq 1/(2*si)
%     wn = 2*wp/si;
%     
%     if (wn > 1)
%         filtdata = data;
%         return
%     end
% 
%     % create the filter coefficients
%     [b,a] = besself(n, wn);
% 
%     % copy
%     filtdata = data;
%     % and replace data portion for all dimensions
%     for i=2:size(filtdata,2)
%         filtdata(:,i) = filtfilt(b,a,data(:,i));
%     end

end
