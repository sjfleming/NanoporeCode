function [ d, h ] = cbfload(filename, range)
    %CBFLOAD Loads data from a CrampEx-generated file in a range
    %   cbfload(filename, range) - Load range=[start,end] of points, 0-based
    %   cbfload(filename, 'info') - Load header data only
    
    d = [];

    % open the file and read header
    fid = fopen(filename,'r');
    if fid == -1
        error(['Could not open ' filename]);
    end
    
    % get total length of file
    d = dir(filename);
    fileSize = d.bytes;

    % now read number of header points, and the header
    nh = uint64(fread(fid,1,'*uint32'));
    hh = fread(fid,nh,'*uint8');
    h = getArrayFromByteStream(hh);
    
    fstart = 4+nh;
    
    % also, how many points and channels?
    h.numChan = numel(h.chNames);
    % this is the total number of points
    h.numTotal = double(fileSize - fstart)/2;
    % this is the number per each channel
    h.numPts = double(h.numTotal/h.numChan);
    
    % load data, or just return header?
    if nargin > 1 && strcmp(range,'info')
        fclose(fid);
        return;
    end
    
    if strcmp(h.type,'IV')
        % if IV-curve, completely ignore range param
        h.numSweeps = numel(h.setVoltages);
        h.numPerSweep = h.numPts/h.numSweeps;
        % load points
        d16 = reshape(fread(fid,h.numTotal,'*int16'),[h.numChan, h.numPerSweep, h.numSweeps]);
        % create output array of doubles
        d = zeros([h.numSweeps, h.numPerSweep, h.numChan]);
        % and copy over channel-by-channel
        for i=1:h.numChan
            d(:,:,i) = double(permute(d16(i,:,:),[3 2 1])) * h.scaleFactors(i);
        end
    else
        % find the right place to seek, and do so
        if nargin<2
            range = [0 h.numPts];
        end
        
        % go to start
        fseek(fid, fstart+range(1)*2*h.numChan, 'bof');
        % how many points and array size
        npts = range(2) - range(1);
        sz = [h.numChan, npts];
        % now read
        d16 = fread(fid, sz, '*int16');
        % create output
        d = zeros(fliplr(size(d16)));
        % and copy/scale
        for i=1:h.numChan
            d(:,i) = double(d16(i,:))*h.scaleFactors(i);
        end
    end

    % all done
    fclose(fid);
end

