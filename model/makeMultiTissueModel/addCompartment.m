function superModel = addCompartment(superModel, compId, compartmentName)
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(superModel,'both');
    identifyMetabolites = sum(abs(superModel.S(:,exchangeRxnsIndexes)),2)>0;
    allMets = superModel.metNames(identifyMetabolites);
    allMets = unique(allMets);
        
    bloodModel.mets =  allMets;
    bloodModel.metNames = allMets;
    bloodModel.rxns =  strcat(bloodModel.mets, '_Exchange');
    bloodModel.rxnNames = bloodModel.rxns;
    conect = [];
    
    %Expand S matrix
    oldSize = size(superModel.S);
    newSize = oldSize + length(allMets);
    newS = sparse(newSize(1), newSize(2));
    
    newS(1:oldSize(1), 1:oldSize(2)) = superModel.S;    

    for i = 1:length(exchangeRxnsIndexes)
        currentMet = find(superModel.S(:, exchangeRxnsIndexes(i)));
        currentMetName = superModel.metNames(currentMet);
        currentMetNr = find(ismember(allMets, currentMetName));
        currentMetRow = currentMetNr + oldSize(1);
        currentMetCol = currentMetNr + oldSize(2);
        newS(currentMetRow, exchangeRxnsIndexes(i)) = 1;
        newS(currentMetRow, currentMetCol) = -1;
    end    
    
    superModel.S = newS;
    bloodModel.name = compartmentName;
    bloodModel.id = compId;
    bloodModel.c = zeros(length(bloodModel.mets),1);
    bloodModel.lb = -1000 * ones(length(bloodModel.rxns),1);
    bloodModel.ub = 1000 * ones(length(bloodModel.rxns),1);
    bloodModel.rev = ones(length(bloodModel.rxns),1);
    bloodModel.b = zeros(length(bloodModel.mets),1);  
    bloodModel.metComps = ones(length(bloodModel.mets),1);
    bloodModel.comps = 's';
    superModel = appendStruct(superModel, bloodModel);
end

