function seq = load_sequence(fn)

    % read sequence from text file
    fid = fopen(fn,'r');
    seqstr = fscanf(fid,'%s');
    fclose(fid);
    % and convert to numerical representation
    seq = nt2int(seqstr)';
end

