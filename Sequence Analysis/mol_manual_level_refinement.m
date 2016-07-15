function mol_manual_level_refinement(mol)
% allows the user to manually refine levels that have been found.
% takes in a molecule object.  modifications are saved in that object.
% Stephen Fleming
% 7/6/16

    % set up a PoreView object to do this
    pv = pv_load_molecule(mol);
    
    % add a 1kHz filtered signal on top
    filtname = 'Low-pass Bessel (1000 Hz)';
    fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_lpb(d,4,1000),filtname), pv.data, 'UniformOutput', false);
    fsigs = fsigs{1};
    pv.setSignalPanel(1,[2,4]);
    
    % initialize and draw levels
    selectedLevels = [0, 0];
    cmap = [0.9290    0.6940    0.1250; ...
        0 1 1];
    redrawLevels();
    
    % display found levels
    function redrawLevels
        % clears everything we have drawn on there
        pv.clearAxes();
        set(pv.getAxes(1),'ColorOrder',cmap);
        % loop through found levels and draw them, selected in bold
        for i = 1:numel(mol.level_means)
            if any(i==selectedLevels)
                pp = plot(pv.getAxes(1), mol.level_timing(i,:)+mol.start_time, mol.level_means(i)*[1,1]*0.001, 'r-', 'LineWidth', 7);
            else
                pp = plot(pv.getAxes(1), mol.level_timing(i,:)+mol.start_time, mol.level_means(i)*[1,1]*0.001, '-', 'LineWidth', 3);
            end
            % set what happens when you click on a level
            set(pp,'ButtonDownFcn',@(~,~) setSel(i));
        end
    end
    
    % merge selected levels
    function mergeLevels
        % merge every level between the two selected
        if sum(selectedLevels~=0) == 2
            % grab all the data
            trange = [mol.level_timing(min(selectedLevels),1), mol.level_timing(max(selectedLevels),2)];
            data = util.doLoadMoleculeData(mol, 10000, 'pointwise', 10000, trange); % load whole thing
            % update the levels
            mol.level_means = [mol.level_means(1:min(selectedLevels)-1); nanmean(data(:,2)); mol.level_means(max(selectedLevels)+1:end)];
            mol.level_medians = [mol.level_medians(1:min(selectedLevels)-1); nanmedian(data(:,2)); mol.level_medians(max(selectedLevels)+1:end)];
            mol.level_stds = [mol.level_stds(1:min(selectedLevels)-1); nanstd(data(:,2)); mol.level_stds(max(selectedLevels)+1:end)];
            mol.level_timing = [mol.level_timing(1:min(selectedLevels)-1,:); trange; mol.level_timing(max(selectedLevels)+1:end,:)];
        end
    end
    
    % split selected levels into two or more levels
    function splitLevel
        % if we have one level selected, split it
        if sum(selectedLevels~=0) == 1
            lev = selectedLevels(selectedLevels~=0);
            % try splitting automatically
            %levels = find_levels_ks(mol, 2000, 5000, 0.07, mol.level_timing(lev,:));
            if false%numel(levels)>1
                mol.level_means = [mol.level_means(1:lev-1); cellfun(@(x) x.current_mean, levels); mol.level_means(lev+1:end)];
                mol.level_medians = [mol.level_medians(1:lev-1); cellfun(@(x) x.current_median, levels); mol.level_medians(lev+1:end)];
                mol.level_stds = [mol.level_stds(1:lev-1); cellfun(@(x) x.current_std, levels); mol.level_stds(lev+1:end)];
                mol.level_timing = [mol.level_timing(1:lev-1,:); ...
                    cell2mat(cellfun(@(x) [x.start_time, x.end_time]+mol.level_timing(lev,1), levels, 'UniformOutput', false)); mol.level_timing(lev+1:end,:)];
            else
                % if no auto-split, let me split manually
                %in = input('No level splittings found... input one manually? (y/n): ','s');
                if true%strcmp(in,'y')
                    %display('Place cursor 1 at the spot of a level transition, then hit any key.')
                    %pause();
                    %trange = pv.getCursors();
                    [t1,~] = ginput(1);
                    %t = trange(1) - mol.start_time; % level timings are based on t=0 at molecule start
                    t = t1 - mol.start_time;
                    % calculate level mean, median, std
                    data = util.doLoadMoleculeData(mol, 10000, 'pointwise', 10000, mol.level_timing(lev,:)); % load whole thing
                    ind = round((t-mol.level_timing(lev,1))/diff(mol.level_timing(lev,:)) * size(data,1));
                    levels{1}.current_mean = nanmean(data(1:ind,2));
                    levels{1}.current_median = nanmedian(data(1:ind,2));
                    levels{1}.current_std = nanstd(data(1:ind,2));
                    levels{1}.start_time = mol.level_timing(lev,1);
                    levels{1}.end_time = t;
                    levels{2}.current_mean = nanmean(data(ind:end,2));
                    levels{2}.current_median = nanmedian(data(ind:end,2));
                    levels{2}.current_std = nanstd(data(ind:end,2));
                    levels{2}.start_time = t;
                    levels{2}.end_time = mol.level_timing(lev,2);
                    if size(levels,1)<size(levels,2)
                        levels = levels';
                    end
                    % update molecule
                    mol.level_means = [mol.level_means(1:lev-1); cellfun(@(x) x.current_mean, levels); mol.level_means(lev+1:end)];
                    mol.level_medians = [mol.level_medians(1:lev-1); cellfun(@(x) x.current_median, levels); mol.level_medians(lev+1:end)];
                    mol.level_stds = [mol.level_stds(1:lev-1); cellfun(@(x) x.current_std, levels); mol.level_stds(lev+1:end)];
                    mol.level_timing = [mol.level_timing(1:lev-1,:); ...
                        cell2mat(cellfun(@(x) [x.start_time, x.end_time], levels, 'UniformOutput', false)); mol.level_timing(lev+1:end,:)];
                end
            end
        end
    end
    
    % set functions for what happens when levels are selected
    function setSel(j)
        if sum(nonzeros(selectedLevels)) == 0 % no selection yet
            selectedLevels(1) = j;
            redrawLevels();
        elseif selectedLevels(1) ~= 0 && j~=selectedLevels(1) % one level selected
            selectedLevels(2) = j;
            redrawLevels();
        end
    end
    
    % set keyboard functions
    function keyFn(e)
        
        switch e.Character
            
            case 'z'
                % save molecule information
                mol.save;
            
            case 'd'
                % deselect levels
                selectedLevels = [0, 0];
                redrawLevels();
                
            case 'r'
                % redraw
                redrawLevels();
            
            case 'm'
                % merge levels
                mergeLevels();
                selectedLevels = [0, 0];
                redrawLevels();
            
            case 's'
                % split levels, if one is selected
                splitLevel();
                selectedLevels = [0, 0];
                redrawLevels();
            
            case 'l'
                % redo level analysis
                % get parameters from user
                str = inputdlg('Enter desired filter frequency, sampling frequency, and p-value:      ','PoreView',1,{'5000, 10000, -30'});
                strs = strsplit(str{1});
                if numel(strs) ~= 3
                    return;
                else
                    filter = str2double(strs{1}(1:end-1));
                    sampling = str2double(strs{2}(1:end-1));
                    p = str2double(strs{3});
                    display('Level finding with parameters:')
                    display(['f = ' num2str(filter) ', s = ' num2str(sampling) ', p = ' num2str(p)])
                end
                mol.do_level_analysis('filter', filter, 'sample', sampling, 'plevels', p);
                redrawLevels();
            
            case 'p'
                % plot in a print-worthy way
                % if cursors, do those
                tr = pv.getCursors();
                if isempty(tr)
                    % otherwise, do the full view
                    tr = pv.getView();
                end
                util.doPlot(pv,tr);
                
            case 'f'
                % ask user what filters they want to add
                str = inputdlg('Enter desired filter and frequency/param:      ','PoreView',1,{'lpb 1000'});
                strs = strsplit(str{1});
                if numel(strs) < 2
                    return
                elseif numel(strs)==3 % band pass
                    params(1) = str2double(strs{2});
                    params(2) = str2double(strs{3});
                    param = 1;
                else
                    param = str2double(strs{2});
                    if isnan(param) || param <= 0
                        return
                    end
                end

                switch strs{1}
                    case 'lp'
                        filtname = sprintf('Low-pass (%d Hz)', param);
                        fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_lp(d,4,param),filtname), pv.data, 'UniformOutput', false);
                        fsigs = fsigs{1};
                    case 'lpb'
                        filtname = sprintf('Low-pass Bessel (%d Hz)', param);
                        fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_lpb(d,4,param),filtname), pv.data, 'UniformOutput', false);
                        fsigs = fsigs{1};
                    case 'hp'
                        filtname = sprintf('High-pass (%d Hz)', param);
                        fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_hp(d,4,param),filtname), pv.data, 'UniformOutput', false);
                        fsigs = fsigs{1};
                    case 'med'
                        filtname = sprintf('Median (%d pts)', param);
                        fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_med(d,param),filtname), pv.data, 'UniformOutput', false);
                        fsigs = fsigs{1};
                    case 'band'
                        filtname = sprintf('Band-pass (%d Hz)', params);
                        fsigs = arrayfun(@(x) x.addVirtualSignal(@(d) filt_band(d,4,params),filtname), pv.data, 'UniformOutput', false);
                        fsigs = fsigs{1};
                    otherwise
                        return
                end            
                % and replace original signals with new ones
                for i=1:numel(pv.psigs)
                    s = pv.psigs(i).sigs;
                    s(s<=pv.data(1).nsigs+1) = s(s<=pv.data(1).nsigs+1) + fsigs(1) - 2;
                    pv.psigs(i).sigs = s;
                end
                pv.refresh();
        end
    end
    
    % and set our all-important keyboard callback
    pv.setKeyboardCallback(@keyFn);

end