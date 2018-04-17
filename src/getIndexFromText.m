function metaboliteNumbers = getIndexFromText(data, search)
    metaboliteNumbers = zeros(length(search),1);
    for i =1:length(metaboliteNumbers)
        currentNumber = findIndex(data, search{i});
        if isempty(currentNumber)
            metaboliteNumbers(i) = -1;
        else
            metaboliteNumbers(i) = findIndex(data, search{i});
        end
    end
end

