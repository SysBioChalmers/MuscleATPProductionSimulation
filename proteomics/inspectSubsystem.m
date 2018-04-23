function inspectSubsystem(model, subsystem, geneMap, tresh)
    allRxns = find(ismember(model.subSystems, subsystem));
    allRxns = model.rxns(allRxns);
    inspectReactions(model, allRxns, geneMap, tresh)
end

