function [i] = findIndex(haystack, nedle)
    i=find(ismember(haystack,nedle));
end
