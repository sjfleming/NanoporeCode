function iterative_scaling_alignment(mol)
%
% Stephen Fleming, 4/5/16
    
    % get the initial predicted levels from oxford
    levs = abs(mol.level_means);
    hicut = 0.6 * abs(mol.open_pore_current);
    lowcut = 0.15 * abs(mol.open_pore_current);
    [model_levels, model_levels_std] = ...
        get_model_levels_oxford(mol.sequence, levs(levs>lowcut & levs<hicut), abs(mol.open_pore_current), abs(mol.voltage), mol.temp);
    
    % save the initial scaling
    mol.predicted_levels = model_levels';
    mol.predicted_levels_stdev = model_levels_std';
    
    % iterate alignment and scaling until convergence is achieved
    stop_criterion = 0.05;
    lsqdist = [];
    dfit = 1;
    while dfit > stop_criterion
        % do a level alignment
        mol.do_level_alignment;
        % calculate least-squares distance per measured level
        logic = ~isnan(mol.level_alignment.model_levels_measured_mean_currents); % these levels are not missing
        lsqdist(end+1) = sum((mol.level_alignment.model_levels_measured_mean_currents(logic) - mol.predicted_levels(logic)).^2) / sum(logic);
        if numel(lsqdist)>1
            dfit = -(lsqdist(end)-lsqdist(end-1));
        else
            figure(3)
            clf(3)
            errorbar(mol.level_alignment.model_levels_measured_mean_currents(logic), ...
                mol.predicted_levels(logic),mol.level_alignment.model_levels_measured_total_duration(logic),'o')
            figure(1)
            clf(1)
            plot(mol.level_alignment.model_levels_measured_mean_currents,'o-')
            hold on
            plot(mol.predicted_levels,'o-')
        end
        % do a least-squares fit
        f = fit(mol.level_alignment.model_levels_measured_mean_currents(logic),mol.predicted_levels(logic), ...
            'poly1','Weights',mol.level_alignment.model_levels_measured_total_duration(logic),'Robust','bisquare');
        % invert to get the scale and offset corrections and update the
        % predicted levels
        mol.predicted_levels = (mol.predicted_levels-f.p2)/f.p1;
    end
    
    plot(mol.predicted_levels,'o-')
    legend('data','initial model','best model')
    
    figure(2)
    clf(2)
    plot(1:numel(lsqdist),lsqdist,'o-')
    xlabel('iterations')
    ylabel('fit distance')
    hold on
    plot(2:numel(lsqdist),abs(diff(lsqdist,1)),'*-')
    
end