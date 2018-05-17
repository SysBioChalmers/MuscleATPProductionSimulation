function atpRates = getATPrates(model, settings)
    maximumSolution = maximizeATP(model, settings);
    maximumSolution = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
    if maximumSolution>(10^-6)
        atpRates = linspace(0, maximumSolution, settings.timeSteps);
    else
        atpRates = [];
    end
end

