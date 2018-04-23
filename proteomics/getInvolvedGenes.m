function allGenes = getInvolvedGenes(model, rxns)
    allGenes = [];
    temp = model.grRules(rxns);
    for i = 1:length(temp)
        currentStr = strrep(temp{i},'(','');
        currentStr = strrep(currentStr,')','');
        content = strsplit(currentStr, ' or ');
        allGenes = [allGenes content];
    end
    allGenes = unique(allGenes);
end

