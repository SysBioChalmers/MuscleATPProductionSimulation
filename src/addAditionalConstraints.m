function [model, allReactions, allBonds] = addAditionalConstraints(model, file)
    additionalConstraints = importdata(file)';
    allReactions = zeros(length(additionalConstraints.data),1);
    allBonds = additionalConstraints.data;
    for i = 1:length(additionalConstraints.data)
        rxn = strsplit(additionalConstraints.textdata{i}, '_back');
        rxnNr = findIndex(model.rxns, rxn{1});
        allReactions(i) = rxnNr;
        if length(rxn)>1
            allBonds(i) = -allBonds(i);
            model.lb(rxnNr) = allBonds(i);
        else
            model.ub(rxnNr) = allBonds(i);
        end
        
    end
end
