function [vO2, fullSolution] = minimizeVO2(model, settings, ATPrate)
    model = setParam(model, 'lb', settings.primaryObjective, ATPrate);
    model.c = zeros(length(model.lb),1);
    currentRxn = getTransport(model, {'O2'}, 'sb', 's');
    model.c(currentRxn) = 1;
    solution = solveLin(model, true);
    fullSolution = solution.x;
    vO2 = fullSolution(currentRxn);
end

