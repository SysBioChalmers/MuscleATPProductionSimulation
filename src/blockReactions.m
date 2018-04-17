function model = blockReactions(model, rxnNames)

for i = 1:length(rxnNames)
   model.lb(findIndex(model.rxns, rxnNames{i})) = 0;    
   model.ub(findIndex(model.rxns, rxnNames{i})) = 0;
end
    
    
end

