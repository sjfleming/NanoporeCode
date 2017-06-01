function [ filtdata ] = filt_conductance( data, voltage_chan )
    % FILT_CONDUCTANCE calculates conductance
    %   Data is passed with colums of [time, current, voltage, ...], and
    %       returns the same array but with the data filtered.
    %   The columns [2:end], except the voltage_chan, get filtered, including virtual signals.

    % copy to output
    filtdata = data;
    % calculate
    for i = 2:size(data,2)
        if i==voltage_chan
            continue;
        end
        filtdata(:,i) = data(:,i)*1000./data(:,voltage_chan) .* double(abs(data(:,voltage_chan))>2);
    end
    % no time column
    filtdata = filtdata(:,2:end);
end

