function newModel = addReversedReactions(model)
    newModel = model;
    SAConstraint = findIndex(model.mets, 'MassConstraint');
    affectedReactions = find(model.S(SAConstraint, :) ~=0);  
    newModel.S = horzcat(model.S, -model.S(:, affectedReactions));  
    newModel.S(SAConstraint,:) = abs(newModel.S(SAConstraint,:));
    newModel.c = [model.c; -model.c(affectedReactions)];
    newModel.rev = [model.rev; zeros(length(affectedReactions),1)];
    newModel.rev(affectedReactions) = 0;
    
    backNames = strcat(model.rxns(affectedReactions), '_back');
    newModel.rxns = [model.rxns; backNames];
    backNames = strcat(model.rxnNames(affectedReactions), '_back');
    newModel.rxnNames = [model.rxnNames; backNames];
%    newModel.rxnComps = [model.rxnComps; model.rxnComps(affectedReactions)];
    newModel.grRules = [model.grRules; model.grRules(affectedReactions)];
    newModel.rxnGeneMat = vertcat(model.rxnGeneMat, model.rxnGeneMat(affectedReactions, :));
    newModel.subSystems = [model.subSystems; model.subSystems(affectedReactions)];
    newModel.eccodes = [model.eccodes; model.eccodes(affectedReactions)];
    newModel.specificActivity = [model.specificActivity; model.specificActivity(affectedReactions)];
    
    
    newRxns = 1:length(affectedReactions);
    newRxns = length(model.rxns) + newRxns;
    
    oldLb = model.lb(affectedReactions);
    oldUb = model.ub(affectedReactions);  
    newModel.lb(newRxns) = -oldUb;
    newModel.ub(newRxns) = -oldLb;

    totalLb = newModel.lb([affectedReactions, newRxns]);
    totalUb = newModel.ub([affectedReactions, newRxns]);
    totalLb(totalLb<0) = 0;
    totalUb(totalUb<0) = 0;
    newModel.lb([affectedReactions, newRxns]) = totalLb;
    newModel.ub([affectedReactions, newRxns]) = totalUb;
end