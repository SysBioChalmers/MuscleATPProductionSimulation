function maximumSolution = maximizeATP(model, settings)
    model = setParam(model, 'obj', settings.primaryObjective, 1);
    if settings.pfba
        solution = solveLinMin(model, true);
    else
        solution = solveLin(model, true);
    end
    if length(solution.x) == 1
        maximumSolution = zeros(size(model.S,2),1);
    else
        maximumSolution = solution.x;
    end
end

