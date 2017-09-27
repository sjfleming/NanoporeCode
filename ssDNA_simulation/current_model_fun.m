function y = current_model_fun(list, p)

% as = p(6:10);
% gs = p(11:15);
p(1:5) = abs(p(1:5))/sum(p(1:5));
y = zeros(size(list));

for i = 1:numel(list)
    
    for j = 1:size(list{1},1)
        switch list{i}(j)
            case 'A'
                y(i) = y(i) + 52.5778 * p(j);% + as(j);
            case 'C'
                y(i) = y(i) + 52.4166 * p(j);
            case 'G'
                y(i) = y(i) + 54.7827 * p(j);% + gs(j);
            case 'T'
                y(i) = y(i) + 36.0750 * p(j);
        end
        
    end
    
end

end