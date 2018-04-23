function  plotAllSubs(model, geneMap, tresh)
    allSubsystems = unique(model.subSystems);
    results = zeros(length(allSubsystems),1);
    for i = 1:length(allSubsystems)
       allRxns = ismember(model.subSystems, allSubsystems{i});
       involvedGenes = getInvolvedGenes(model, allRxns);
       results(i) = sumOfGenes(involvedGenes);
    end
    [results, indx] = sort(results, 'descend');
    allSubsystems = allSubsystems(indx);
    allSubsystems(results<tresh) = [];
    results(results<tresh) = [];
    figure()
    plotPie(results, allSubsystems)
    allSubsystems
end

