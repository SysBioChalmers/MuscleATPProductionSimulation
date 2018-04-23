function printGeneList(geneMap, geneList)

    results = zeros(length(geneList),1);
    
    for i = 1:length(geneList)
        results(i) = geneMap(geneList{i});
    end
    
    [results indx] = sort(results, 'descend');
    geneList = geneList(indx);
    
    for i = 1:length(geneList)
        fprintf('%s\t%2.5f\n', geneList{i}, results(i));
    end    

end

