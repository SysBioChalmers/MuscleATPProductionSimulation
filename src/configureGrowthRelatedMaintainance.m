function model = configureGrowthRelatedMaintainance(model, ATP, objective, growthRelatedMaintainance)
    reaction = findIndex(model.rxns, objective);
    metaboliteNames = modifyMetNames(model);
    metabolite = findIndex(metaboliteNames, growthRelatedMaintainance);
    model.S(metabolite, reaction) = -ATP;
end

