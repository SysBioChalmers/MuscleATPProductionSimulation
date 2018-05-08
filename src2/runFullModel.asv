function [fullSolution, growthRates] = runFullModel(model, settings)


inMedia = settings.inMedia;
constraints = settings.inValues;
outMedia = settings.outMedia;
outConstraints = settings.outConstraints;
nrTimePoints = settings.timeSteps;
objectiveFunction = settings.primaryObjective;
settings.minimizationRxns = [outMedia(2) inMedia(3)];
settings.minimizationCoeficients = [-1 0.5];



[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');


%Remove ProtPool

model.foodRxns = getBounds(model, inMedia);
reactionIn= model.foodRxns(model.foodRxns ~= -1);

model.foodRxns = getBounds(model, outMedia);
reactionOut = model.foodRxns(model.foodRxns ~= -1);


model = setParam(model, 'lb', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', exchangeRxnsIndexes, 0);

model = setParam(model, 'lb', reactionIn, constraints);
model = setParam(model, 'ub', reactionOut, outConstraints);

model = setParam(model, 'obj', objectiveFunction, 1);
solution = solveLin(model, true);
maximumSolution = abs(solution.f);

if maximumSolution>(10^-6)
    growthRates = linspace(0.01, maximumSolution, nrTimePoints);


    model.c = zeros(length(model.lb),1);
    for i = 1:length(settings.minRxns)
        currentRxn = findIndex(model.rxns, settings.minRxns{i});
        model.c(currentRxn) = settings.minVal(i);
    end
    
    if isfield(settings, 'QminRxns')
        quadratic = true;
        model.qc = zeros(length(model.c),1);
        for i = 1:length(settings.QminRxns)
            currentRxn = findIndex(model.rxns, settings.QminRxns{i});
            model.qc(currentRxn) = settings.QminVal(i);
        end
    else
        quadratic = false;
    end
    
    
    %Minimize Oxygen
    %model = setParam(model, 'obj', reactionIn(3), 1);    

    %Minimize Oxygen and Glycogen
    %model = setParam(model, 'obj', reactionIn([1, 3]), [1 0.5]);
    
    %Minimize CO2
    %model = setParam(model, 'obj', reactionOut(2), -1);

    fullSolution =  zeros(nrTimePoints, length(model.rxns));



    for i = 1:nrTimePoints
        %A too strikt bond on ATP production rate may give false
        %infeasibility
        model = setParam(model, 'lb', objectiveFunction, growthRates(i));
        model = setParam(model, 'ub', objectiveFunction, growthRates(i)*1.001);
        
        
        %Constrain blood floow as well
        %model.ub(findIndex(model.rxns, 'BloodCost')) = growthRates(i)*0.155;

        if quadratic == true
            solution = quadOpt(model);
        else
            solution = solveLin(model, true);
            %solution = solveLinMin(model, true);
        end

        %solution = minimizeTransport(model, solution, 'sm', 'sb');

        %solution = solveLP(model);

        if length(solution.x)>1
            fullSolution(i,:) = solution.x;
        end
    end
else
    fullSolution=0;
    growthRates=0;
end
    %printFluxesAvlant(model, solution.x, true)

end


function newSolution = minimizeTransport(model, solution, from, to)
    mTransp = find(getTransportReactions(model, from, to));
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
    model.lb(exchangeRxnsIndexes) = solution.x(exchangeRxnsIndexes);
    model.ub(exchangeRxnsIndexes) = solution.x(exchangeRxnsIndexes);
    emptyTransport =  abs(solution.x(mTransp))<10^-5;
    mTransp(emptyTransport) = [];
    objective = solution.x(mTransp);
    objective = sign(objective);
    model = setParam(model, 'obj', mTransp, objective);
    newSolution = solveLin(model, true);
end
