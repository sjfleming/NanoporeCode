function hps = betterbar(xs, ns, spc)

    nbars = size(ns,2);
    % total distance between bars
    xspc = xs(2)-xs(1);
    % spc is extra space between sets of bars, fractional
    if nbars > 1
        delt = 0.5*xspc*(1-spc)*(1-1/nbars);
        delts = linspace(-delt,delt,nbars);
    else
        delts = 0;
    end
    
    
    colors = [194 224 250; 210 226 139; 255 154 200]/255.0;

    for j=1:nbars
        bar(xs+delts(j),ns(:,j),(1-spc)/nbars,'FaceColor',colors(j,:))
        hold on
    end
    
end

