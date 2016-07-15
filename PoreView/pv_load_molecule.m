function pv = pv_load_molecule(mol)
% Load the files corresponding to a molecule object into a PoreView session

% Stephen Fleming
% 7/7/16

    % get file info from molecule
    firstfile = mol.start_file;
    lastfile = mol.end_file;
    
    % create a PoreView and load the first file
    pv = PoreView(firstfile);
    
    % figure out all the files involved
    files = {firstfile};
    if ~strcmp(firstfile,lastfile)
        firstfilenum = str2double(firstfile(end-7:end-4));
        lastfilenum = str2double(lastfile(end-7:end-4));
        prefix = firstfile(1:end-8);
        suffix = firstfile(end-3:end);
        if lastfilenum-firstfilenum+1 > 5
            in = input(['Trying to load ' num2str(lastfilenum-firstfilenum+1) ' files.  Proceed? (y/n): '],'s');
            if strcmp(in,'n')
                pv = [];
                return;
            end
        end
        for i = (firstfilenum+1):lastfilenum
            files = cat(1, files, [prefix sprintf('%04d',i) suffix]);
        end
        
        % load all these additional files
        for i = 2:numel(files)
            pv.catFile(files{i});
        end
        
    end
    
    % zoom into a range around the molecule and set cursors
    trange = mol.start_time + [0, mol.level_timing(end,2)];
    pv.setCursors(trange);
    dt = diff(trange)/50;
    pv.setView([trange(1)-dt trange(2)+dt])

end