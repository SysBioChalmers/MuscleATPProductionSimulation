function [param, paramComments] = loadParm(subject, print)
    data = importdata(['sampleData/' subject '/param.txt'], '\t');
    param = [];
    for i = 1:length(data)
       curStr = strsplit(data{i},'\t');
       curVal = str2num(curStr{2});
       param.(curStr{1}) = curVal;
       paramComments.(curStr{1}) = curStr{3};
       if print
          fprintf('%s\t%2.2f\t%s\n', curStr{1}, curVal, curStr{3}) 
       end
    end
end

