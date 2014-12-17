function [ d, h ] = fast5load(filename, range, channels)
    %FAST5LOAD Loads raw data from a MinKNOW-generated fast5 file
    %   fast5load(filename, range, channels) - Load range=[start,end] of
    %       points, in the specified channels (1-512)
    %   fast5load(filename, 'info') - Load header data only, return which
    %       channels are recorded in the file


    d = [];
    
    % unlike other ones, we're only returning header here if info requested
    % otherwise, we assume channels and range is valid...!
    if strcmp(range,'info')
        % header request
        h = [];
        % check the channels
        numpts = 0;
        digitization = 0;
        mrange = 0;
        chans = [];
        h.offset = zeros(512,1);
        for chan = 1:512
            chanstr = ['/Raw/Channel_', num2str(chan), '/Signal'];
            chanmeta = ['/Raw/Channel_', num2str(chan), '/Meta'];
            try
                s = h5info(filename, chanstr);
            catch
                % this channel isn't present, so skip it
                continue
            end
            % channel exists
            chans(end+1) = chan;
            numpts = s.Dataspace.Size;
            h.offset(chan) = h5readatt(filename, chanmeta, 'offset');
            digitization = h5readatt(filename, chanmeta, 'digitisation');
            mrange = h5readatt(filename, chanmeta, 'range');
        end
        
        h.numPts = numpts;
        % min and max channels in the original file
        h.minChan = min(chans);
        h.maxChan = max(chans);
        
        % digitization and range
        h.digitization = digitization;
        h.range = mrange;
        
        samplerate = h5readatt(filename, chanmeta, 'sample_rate');
        
        h.si = 1.0/samplerate;
        
        return;
    end

    % no header, just return data
    d = zeros(range(2)-range(1),numel(channels));
    
    for i=1:numel(channels)
        chan = channels(i);
        chanstr = ['/Raw/Channel_', num2str(chan), '/Signal'];
        chanmeta = ['/Raw/Channel_', num2str(chan), '/Meta'];
        
        % get various attributes for each channel
        % yes, this slows it down a bit, but oh well
        % (abf load does it too!)
        offset = h5readatt(filename, chanmeta, 'offset');
        digitization = h5readatt(filename, chanmeta, 'digitisation');
        adcrange = h5readatt(filename, chanmeta, 'range');

        % read corresponding channel's points
        d(:,i) = h5read(filename, chanstr, range(1)+1, range(2)-range(1));
        % and scale it according to metadata
        d(:,i) = (d(:,i)+offset)*adcrange/digitization;
    end
end

