function mol = batch_sliding_analysis(date,filenums,temp,KCl_molarity,ADPNP_molarity)

% Stephen Fleming
% 20150210

% analysis of sliding data in automatic batches

% initialize
addpath('/Users/Stephen/GitHub/NanoporeCode/PoreView/');
folder = '/Users/Stephen/Documents/Stephen/Research/Data/Biopore/';

% example:
% date = '20160121';
% filenums = 7:12;

files = cell(0);
mol = cell(0);
unfinished_event_in_last_file = false;

for k = 1:numel(filenums)
    open_pore = nan;
    files{k} = [date '/' date(1:4) '_' date(5:6) '_' date(7:8) '_' sprintf('%04d',filenums(k)) '.abf'];
    if (ishandle(1))
        close(1)
    end
    pv = pv_launch([folder files{k}]);
    sigdata = pv.data;
    %sigdata = SignalData([folder files{k}]);
    display([files{k} ' : '])
    
    % find the good events
    events = util.doFindEvents(sigdata,1);
    goodEvent = false(1,numel(events));
    for j = 1:numel(events)
        events{j}.start_file = [folder files{k}];
        events{j}.end_file = [folder files{k}];
        events{j}.temp = temp;
        events{j}.ADPNP_molarity = ADPNP_molarity;
        events{j}.KCl_molarity = KCl_molarity;
        
        % allow me to manually add and reject events
        pv.setCursors([events{j}.start_time, events{j}.end_time]);
        pause();
        tr = pv.getCursors();
        if ~(tr(1)==events{j}.start_time && tr(2)==events{j}.end_time)
            if strcmp(input('Good event? (y/n): ','s'),'y')
                events{j}.start_time = tr(1);
                events{j}.end_time = tr(2);
                events{j}.start_ind = tr(1)/pv.data.si;
                events{j}.end_ind = tr(2)/pv.data.si;
                events{j}.ended_manually = strcmp(input('Is this event ended manually? (y/n): ','s'),'y');
                events{j}.continues_past_end_of_file = strcmp(input('Does this event continue past the end of the file? (y/n): ','s'),'y');
                goodEvent(j) = true;
                if isnan(open_pore) || isempty(open_pore)
                    display('Set cursors on a segment of open pore current.')
                    pause();
                    open_pore = mean(pv.data.get(pv.getCursors()/pv.data.si,2))*1000;
                    events{j}.open_pore_current = open_pore;
                    events{j}.voltage = mean(pv.data.get(pv.getCursors()/pv.data.si,3));
                end
            else
                goodEvent(j) = false;
            end
        else
            goodEvent(j) = true;
        end
    end
    
    % missed any events?
    display('That was the last auto-found event.')
    while strcmp(input('Input another event manually? (y/n): ','s'),'y')
        j = j+1;
        if isempty(j)
            j = 1;
        end
        display('Set cursors to beginning and end of event, then hit any key.')
        pause();
        tr = pv.getCursors();
        if strcmp(input('Good event? (y/n): ','s'),'y')
            events{j}.start_time = tr(1);
            events{j}.end_time = tr(2);
            events{j}.start_ind = tr(1)/pv.data.si;
            events{j}.end_ind = tr(2)/pv.data.si;
            events{j}.ended_manually = strcmp(input('Is this event ended manually? (y/n): ','s'),'y');
            events{j}.continues_past_end_of_file = strcmp(input('Does this event continue past the end of the file? (y/n): ','s'),'y');
            events{j}.start_file = [folder files{k}];
            events{j}.end_file = [folder files{k}];
            events{j}.temp = temp;
            events{j}.ADPNP_molarity = ADPNP_molarity;
            events{j}.KCl_molarity = KCl_molarity;
            goodEvent(j) = true;
            display('Set cursors on a segment of open pore current.')
            pause();
            events{j}.open_pore_current = mean(pv.data.get(pv.getCursors()/pv.data.si,2))*1000;
            events{j}.voltage = mean(pv.data.get(pv.getCursors()/pv.data.si,3));
        else
            goodEvent(j) = false;
        end
    end
    
    events = events(goodEvent);
    for i = 1:numel(events)
        % combine first with last from previous file, if necessary
        if (numel(mol)>=1 && i == 1 && ...
                unfinished_event_in_last_file ...
                && events{i}.start_ind < 100)
            mol{end}.addData(events{i});
            display(['Last molecule in ' mol{end}.start_file(end-18:end)])
            display(['combined with first molecule in ' events{i}.start_file(end-18:end)])
            mol{end}.save;
        else
            mol{end+1} = molecule(events{i});
            mol{end}.save;
        end
    end
    % check if last molecule is unfinished in this file
    if ~isempty(i) && events{i}.continues_past_end_of_file;
        unfinished_event_in_last_file = true;
    else
        unfinished_event_in_last_file = false;
    end
    %events = [events, events_temp(goodEvent)];
    %events = [events, events_temp];
    clear events_temp sigdata j pv;
    
end

% check for incomplete molecules and combine with recording in the next file
i = 1;
a = 1;
while i <= numel(events)
    mol{a} = molecule(events{i});
    if (events{i}.continues_past_end_of_file) && numel(events)>i
        if (str2num(events{i+1}.start_file(end-6:end-4)) == str2num(events{i}.end_file(end-6:end-4))+1 ...
                && events{i+1}.start_ind < 100)
            mol{a}.addData(events{i+1});
            display(['Last molecule in ' mol{a}.start_file(end-18:end)])
            display(['combined with first molecule in ' events{i+1}.start_file(end-18:end)])
            %display(mol{i})
            i = i+2;
        else % molecule didn't end in first file, but it doesn't continue in next either, it's just cut off...
            i = i+1;
        end
    else
        i = i+1;
    end
    mol{a}.save;
    a = a+1;
end

% % do the level analysis and save data
% for i = 1:numel(mol)
%     filter = 1000;
%     sampling = 5000;
%     p = -15;
%     display(['Molecule ' num2str(i) ' level analysis']);
%     mol{i}.do_level_analysis('sample', sampling, 'plevels', p, 'filter', filter);
%     mol{i}.save;
% end
% 
% % add in the sequence data and align levels to predictions
% 
% % adapterSF with the hairpin
% % 'R' = abasic
% seq = 'RRRRRTTTTTTTTTTTTGGGAAATTTTTGGGAAATTTTCGATCACTGGAACTTTACAAGGAATTTCCTGTGAAGCTGCCGAGGTTTGACGCGARRRACATGACGGGATGCGGAATCTTTTGATTCCGCATCCCGTCATGTTGCTCGCGTCAAACCTCGGCAGCTTCACAGGAAATTCCTTGTAAAGTTCCAGTGATCGAAAATTTCCCAAAAATTTCCCTTTGAGGCGAGCGGTCAA';
% %seq = 'RRRRRTTTTTTTTTTTTGGGAAATTTTTGGGAAATTTTCGATCACTGGAACTTTACAAGGAATTTCCT';
% 
% for i = 1:numel(mol)
%     display(['Molecule ' num2str(i) ' level alignment']);
%     
%     try
%         % get the predicted levels from oxford
%         levs = abs(mol{i}.level_means);
%         hicut = 0.6 * abs(mol{i}.open_pore_current);
%         lowcut = 0.15 * abs(mol{i}.open_pore_current);
%         %     [model_levels, model_levels_std] = get_model_levels_oxford(seq, levs(levs>lowcut & levs<hicut), abs(mol{i}.open_pore_current), abs(mol{i}.voltage), mol{i}.temp);
%         [model_levels, model_levels_std] = get_model_levels_my_M2(seq, levs(levs>lowcut & levs<hicut));
%         mol{i}.predicted_levels = model_levels';
%         mol{i}.predicted_levels_stdev = model_levels_std';
%         mol{i}.sequence = seq;
%         
%         % do the alignment of measured levels to predictions
%         mol{i}.do_level_alignment;
%         
%         % save everything
%         mol{i}.save;
%     catch ex
%         display(['Skipping alignment of molecule ' num2str(i) ': problems.']);
%     end
% end

display('Complete.')

end

