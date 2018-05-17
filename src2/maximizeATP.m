function maximumSolution = maximizeATP(model, settings)
    model = setParam(model, 'obj', settings.primaryObjective, 1);
    if settings.pfba
        solution = solveLinMin(model, true);
    else
        solution = solveLin(model, true);
    end
    maximumSolution = solution.x;
end

