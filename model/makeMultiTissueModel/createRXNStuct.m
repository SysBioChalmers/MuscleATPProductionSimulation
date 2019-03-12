function rxnsToAdd = createRXNStuct(model, name, equation, lb, ub, subsystem)
    rxnsToAdd.rxns{1} = name;
    rxnsToAdd.rxnNames{1} = '';
    rxnsToAdd.equations{1} = equation;
    rxnsToAdd.lb = lb;
    rxnsToAdd.ub = ub;
    rxnsToAdd.subSystems{1} = subsystem;
end