function model = addUncoupling(model)
rxnStruct.rxns{1} = 'Uncoupling';
rxnStruct.equations{1} = 'H+[c] => H+[m]';
rxnStruct.rxnNames{1} ='Uncoupling';
rxnStruct.lb(1) = 0;
rxnStruct.ub(1) = 1000;
rxnStruct.subSystems{1} = 'Oxidative phosphorylation';
rxnStruct.rxnComps(1) = 4;

model = addRxns(model, rxnStruct, 3);
end

