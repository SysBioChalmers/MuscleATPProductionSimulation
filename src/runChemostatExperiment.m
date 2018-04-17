function fullSolution = runChemostatExperiment(model, growthRates, objectiveFunction, minimizeRxn)
fluxThresh = 10^-3;

massConstraintRow = findIndex(model.mets, 'MassConstraint');

fullSolution =  zeros(length(growthRates), length(model.rxns));
weightEstimate = zeros(length(growthRates), 1);
weightRow = full(model.S(massConstraintRow,:)); 



for i = 1:length(growthRates)
    model = setParam(model, 'lb', objectiveFunction, growthRates(i));
    model = setParam(model, 'ub', objectiveFunction, growthRates(i));

    %solution = solveLP(model, 1);
    solution = solveLinMin(model, true);
    if length(solution.x)>1
        fullSolution(i,:) = solution.x;

        weightEstimate(i) = weightRow * solution.x;
    end
    %printFluxesAvlant(model, solution.x, false)
end
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);

solution = solveLin(model)


end

