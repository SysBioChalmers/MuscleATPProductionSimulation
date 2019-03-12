function atpRates = getATPrates(model, settings)
    maxSol = maximizeATP(model, settings);
    maxSol = maxSol(findIndex(model.rxns, settings.primaryObjective));
    if maxSol>(10^-6)
        atpRates = linspace(0, maxSol, settings.timeSteps);
    else
        atpRates = [];
    end
end

