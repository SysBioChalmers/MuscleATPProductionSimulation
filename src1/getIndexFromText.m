function metaboliteNumbers = getIndexFromText(data, search)
% getIndexFromText
% returns the index in data for multiple search strings in  search
%
%   data                a cell array with strings
%   search              a cell array with strings
%   metaboliteNumbers   a list of indexes of the strings in search
%
%   Avlant Nilsson, 2016-05-16
%
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