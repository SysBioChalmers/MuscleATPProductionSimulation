function saveParam(subject, param)
    [oldParam, paramComments] = loadParm(subject, false);
    
    fns = fieldnames(param);
    
    fileID = fopen(['sampleData/' subject '/param.txt'],'w');
    
    
    for i = 1:length(fns)
       curValue = param.(fns{i});
       curComment = paramComments.(fns{i});
       fprintf(fileID,'%s\t%2.2f\t%s\n', fns{i}, curValue, curComment);
    end
    fclose(fileID);
end