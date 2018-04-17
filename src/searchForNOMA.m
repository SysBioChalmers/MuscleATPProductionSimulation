function fullSolution = searchForNOMA(model, growthRates, fullSolution)
threshold = 10^-3;
massConstraintRow = findIndex(model.mets, 'MassConstraint');
objectiveFunction = {'human_ATPMaintainance'};
weightValue = model.b(massConstraintRow,2);
weightRow = full(model.S(massConstraintRow,:));

maxFluxes = max(fullSolution);

weightEstimate = sum(maxFluxes .* weightRow);

weightFactor = weightEstimate/weightValue;
maxFluxes = maxFluxes./weightFactor;

reactionsWithMass = find(weightRow>threshold);

model.ub(reactionsWithMass) = maxFluxes(reactionsWithMass);

model.b(massConstraintRow,:) = [0 1000];

[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
fullSolution =  zeros(length(growthRates), length(model.rxns));

for i = 1:length(growthRates)
    model = setParam(model, 'lb', objectiveFunction, growthRates(i));
    model = setParam(model, 'ub', objectiveFunction, growthRates(i));

    solution = solveLinMin(model, true);
    if length(solution.x)>1
        results(i,:) = solution.x(exchangeRxnsIndexes);
        fullSolution(i,:) = solution.x;
    end
    sum(solution.x .* weightRow')
end

end

