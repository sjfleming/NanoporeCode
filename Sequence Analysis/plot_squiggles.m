function plot_squiggles( levels, durations, t0, axes )
    % plots squiggle data
    
    if nargin < 4
        pv = PoreView();
        axes = pv.getAxes(1);
    end
    if nargin < 3
        t0 = 0;
    end
    if nargin < 2
        durations = zeros(size(levels)) + 0.1;
    end
    
    if numel(durations) ~= numel(levels)
        error('Levels and times mismatch!');
    end
    
    % Build a new color map
    CO = [  0.0 0.2 0.7;...
            0.0 0.5 0.0;...
            0.7 0.1 0.1;...
            0.0 0.6 0.6;...
            0.5 0.0 0.5;...
            0.6 0.6 0.3  ];
    
    for i=1:size(levels,2)
        ts = t0 + [0; cumsum(durations(:,i))];
        l2 = doublemat(levels(:,i));
        t2 = doublemat(ts);
    
        plot(axes, t2(2:end-1), l2, 'Color', CO(i,:), 'LineWidth', 1.5);
    end
    
    if nargin < 4
        pv.psigs(1).setY([min(rshape(levels)),max(rshape(levels))]);
        pv.setView([min(t2) max(t2)]);
    end
end

