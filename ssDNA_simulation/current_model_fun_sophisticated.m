function y = current_model_fun_sophisticated(list, p)

As = p(5+(1:31));
Gs = p(5+31+(1:31));
weights = abs(p(1:5))/sum(p(1:5)); % the simple averaging weights

allcombosA = regexprep(regexprep(string(dec2bin(1:31)),'1','A'),'0','*');
allcombosG = regexprep(regexprep(string(dec2bin(1:31)),'1','G'),'0','*');
y = zeros(size(list));

for i = 1:numel(list)
    
    % the simple weighted average level
    for j = 1:numel(list{i})
        switch list{i}(j)
            case 'A'
                y(i) = y(i) + 52.5778 * weights(j);
            case 'C'
                y(i) = y(i) + 52.4166 * weights(j);
            case 'G'
                y(i) = y(i) + 54.7827 * weights(j);
            case 'T'
                y(i) = y(i) + 36.0750 * weights(j);
        end
    end
    
    % the corrections for As
    strA = regexprep(list{i}','[^A]','*'); % keeps only the As, rest are *
    if ~strcmp(strA,'*****')
        y(i) = y(i) + As(find(strcmpi(strA,allcombosA),1,'first'));
    end
    
    % the corrections for Gs
    strG = regexprep(list{i}','[^G]','*'); % keeps only the Gs, rest are *
    if ~strcmp(strG,'*****')
        y(i) = y(i) + Gs(find(strcmpi(strG,allcombosG),1,'first'));
    end
    
end

end