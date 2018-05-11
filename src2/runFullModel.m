function [atpRates, fullSolution] = runFullModel(model, settings)

inMedia = settings.inMedia;
constraints = settings.inValues;
outMedia = settings.outMedia;
outConstraints = settings.outConstraints;
nrTimePoints = settings.timeSteps;
objectiveFunction = settings.primaryObjective;
settings.minimizationRxns = [outMedia(2) inMedia(3)];
settings.minimizationCoeficients = [-1 0.5];
parsimonious = settings.pfba;


[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');


%Remove ProtPool
reactionIn= getBounds(model, inMedia);
reactionOut = getBounds(model, outMedia);


model = setParam(model, 'lb', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', exchangeRxnsIndexes, 0);

model = setParam(model, 'lb', reactionIn, constraints);
model = setParam(model, 'ub', reactionOut, outConstraints);

model = setParam(model, 'obj', objectiveFunction, 1);
solution = solveLin(model, true);
maximumSolution = abs(solution.f);

if maximumSolution>(10^-6)
    atpRates = linspace(0, maximumSolution, nrTimePoints);


    model.c = zeros(length(model.lb),1);
    for i = 1:length(settings.minRxns)
        currentRxn = findIndex(model.rxns, settings.minRxns{i});
        model.c(currentRxn) = settings.minVal(i);
    end
        
    fullSolution =  zeros(nrTimePoints, length(model.rxns));

    for i = 1:nrTimePoints
        %A too strikt bond on ATP production rate may give false
        %infeasibility
        model = setParam(model, 'lb', objectiveFunction, atpRates(i));
        model = setParam(model, 'ub', objectiveFunction, atpRates(i)*1.001);
        
        if parsimonious == true
            solution = solveLinMin(model, true);
        else
            solution = solveLin(model, true);
        end

        if length(solution.x)>1
            fullSolution(i,:) = solution.x;
        end
    end
else
    fullSolution=0;
    atpRates=0;
end

end
