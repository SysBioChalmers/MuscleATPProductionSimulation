function model = setObjectiveFunction(model, settings)
    model.c = zeros(length(model.lb),1);
    for i = 1:length(settings.minRxns)
        currentRxn = findIndex(model.rxns, settings.minRxns{i});
        model.c(currentRxn) = settings.minVal(i);
    end
end

