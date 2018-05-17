function [i] = findIndex(haystack, nedle)
% findIndex
% Searches for a string (nedle) in a string cell array (haystack)
%
%   haystack   a cell array with strings 
%   nedle      a string
%   [i]        the index in the cell array which matches the string
%
%   Avlant Nilsson, 2016-05-16
%
    i=find(ismember(haystack,nedle));
end
