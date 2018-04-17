function model = removeSubsystems(model, subsystems)
    for i = 1:length(subsystems)
        model.lb(ismember(model.subSystems, subsystems{i})) = 0;
        model.ub(ismember(model.subSystems, subsystems{i})) = 0;
    end
end

