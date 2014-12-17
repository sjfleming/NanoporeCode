% Stephen Fleming
% 4/14/14

function filterData = discretize(data,D,cuts,n)

% Discretizes dataset based on inputs:
% D = array of acceptable levels
% cuts = array of cutoffs that determine limits of the levels
% n = number of elements to average over to determine discrete level
% note: numel(cuts) = numel(D)-1    is a requirement

filterData = zeros(numel(data),1);
i = 1;
n = 1000;
while i+jump < numel(data)
    
    if m <= 25
        level = 19;
    elseif (m>25 && m<=37)
        level = 30;
    elseif (m>37 && m<=51)
        level = 41;
    elseif m>51
        level = 55;
    end
    
    filterData(i:i+n,1) = level;
    
    if mod(round(i/numel(data)*10000)/100,10)==0
        display([num2str(round(i/numel(data)*100)) '% complete'])
    end
    
    i = i+n;
    
end

end
