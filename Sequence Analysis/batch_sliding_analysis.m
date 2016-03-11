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

events = cell(0);
files = cell(0);

for k = 1:numel(filenums)
    files{k} = [date '/' date(1:4) '_' date(5:6) '_' date(7:8) '_' sprintf('%04d',filenums(k)) '.abf'];
    if (ishandle(1))
        close(1)
    end
    pv = PoreView([folder files{k}]);
    sigdata = pv.data;
    display([files{k} ' : '])
    
    % find the good events
    events_temp = util.doFindEvents(sigdata,1);
    goodEvent = false(1,numel(events_temp));
    for j = 1:numel(events_temp)
        events_temp{j}.start_file = [folder files{k}];
        events_temp{j}.end_file = [folder files{k}];
        events_temp{j}.temp = temp;
        events_temp{j}.ADPNP_molarity = ADPNP_molarity;
        events_temp{j}.KCl_molarity = KCl_molarity;
        
        % allow me to manually add and reject events
        pv.setCursors([events_temp{j}.start_time, events_temp{j}.end_time]);
        pause();
        tr = pv.getCursors();
        if ~(tr(1)==events_temp{j}.start_time && tr(2)==events_temp{j}.end_time)
            if strcmp(input('Good event? (y/n): ','s'),'y')
                events_temp{j}.start_time = tr(1);
                events_temp{j}.end_time = tr(2);
                events_temp{j}.start_ind = tr(1)/pv.data.si;
                events_temp{j}.end_ind = tr(2)/pv.data.si;
                events_temp{j}.ended_manually = strcmp(input('Is this event ended manually? (y/n): ','s'),'y');
                events_temp{j}.continues_past_end_of_file = strcmp(input('Does this event continue past the end of the file? (y/n): ','s'),'y');
                goodEvent(j) = true;
            else
                goodEvent(j) = false;
            end
        else
            goodEvent(j) = true;
        end
    end
    
    events = [events, events_temp(goodEvent)];
    clear events_temp sigdata j pv;
    
end

% check for incomplete molecules and combine with recording in the next file
i = 1;
a = 1;
while i <= numel(events)
    mol{a} = molecule(events{i});
    if (events{i}.continues_past_end_of_file)
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

% do the level analysis and save data
for i = 1:numel(mol)
    filter = 2000;
    sampling = 5000;
    p = -25;
    display(['Molecule ' num2str(i) ' level analysis']);
    mol{i}.do_level_analysis(filter, sampling, p);
    mol{i}.save;
end

% add in the sequence data and align levels to predictions

% adapterSF with the hairpin
% 'R' = abasic
%seq = 'RRRRRTTTTTTTTTTTTGGGAAATTTTTGGGAAATTTTCGATCACTGGAACTTTACAAGGAATTTCCTGTGAAGCTGCCGAGGTTTGACGCGARRRACATGACGGGATGCGGAATCTTTTGATTCCGCATCCCGTCATGTTGCTCGCGTCAAACCTCGGCAGCTTCACAGGAAATTCCTTGTAAAGTTCCAGTGATCGAAAATTTCCCAAAAATTTCCCTTTGAGGCGAGCGGTCAA';
seq = 'RRRRRTTTTTTTTTTTTGGGAAATTTTTGGGAAATTTTCGATCACTGGAACTTTACAAGGAATTTCCT';

for i = 1:numel(mol)
    display(['Molecule ' num2str(i) ' level alignment']);
    
    % get the predicted levels from oxford
    levs = abs(mol{i}.level_means);
    hicut = 0.6 * abs(mol{i}.open_pore_current);
    lowcut = 0.1 * abs(mol{i}.open_pore_current);
    model_levels = get_model_levels_oxford(seq, levs(levs>lowcut & levs<hicut), abs(mol{i}.open_pore_current), abs(mol{i}.voltage), mol{i}.temp);
    mol{i}.predicted_levels = model_levels';
    mol{i}.sequence = seq;
    
    % do the alignment of measured levels to predictions
    mol{i}.do_level_alignment;
    
    % save everything
    mol{i}.save;
end

end

