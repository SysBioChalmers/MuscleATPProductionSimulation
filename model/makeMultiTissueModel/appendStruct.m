function superModel = appendStruct(superModel, curModel)
    superModel.rxns = [superModel.rxns; strcat(curModel.rxns, ['_' curModel.id])];
    superModel.mets = [superModel.mets; strcat(curModel.mets, curModel.id)];
    superModel.rxnNames = [superModel.rxnNames; curModel.rxnNames];
    superModel.metNames = [superModel.metNames; curModel.metNames];        
    %superModel.rxnComps = [superModel.rxnComps; curModel.rxnComps];
    superModel.metComps = [superModel.metComps; length(superModel.comps) + curModel.metComps];
    superModel.comps = [superModel.comps; strcat(curModel.comps, curModel.id)];
    superModel.lb = [superModel.lb; curModel.lb];
    superModel.ub = [superModel.ub; curModel.ub];
    superModel.rev = [superModel.rev; curModel.rev];
    superModel.b = [superModel.b; curModel.b];
    superModel.c = [superModel.c; curModel.c];
    superModel.subSystems = [superModel.subSystems; strcat(curModel.subSystems, ['_' curModel.name])];
end