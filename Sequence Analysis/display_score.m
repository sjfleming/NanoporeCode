function widget=display_score()
    % returns a struct/fake class to handle serial output
    % closures for the win!
    
    % first close any open ports
    ser = instrfindall;
    if ~isempty(ser)
        fclose(ser);
    end
    % then try to initialize one
    ser = [];
    try
        ser = serial('COM4','BaudRate',9600);
        fopen(ser);
        pause(2);
    catch
        fprintf(2,'Error connecting to serial!\n');
    end
    % and a function to display the score
    function serialprint(score)
        % output string to serial
        try
            serstr = sprintf('%4.1f',score);
            for i=1:numel(serstr)
                fwrite(ser,serstr(i),'uint8');
            end
        end
    end

    widget = [];
    widget.print = @serialprint;
    widget.cleanup = onCleanup(@() fclose(ser));
end