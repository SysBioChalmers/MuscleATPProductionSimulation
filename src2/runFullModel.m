function [atpRates, fullSolution] = runFullModel(model, settings)
objectiveFunction = settings.primaryObjective;

atpRates = getATPrates(model, settings);

if not(isempty(atpRates))

    model = setObjectiveFunction(model, settings);
        
    fullSolution =  zeros(length(atpRates), length(model.rxns));

    for i = 1:length(atpRates)
        %A too strikt bond on ATP production rate may give false
        %infeasibility
        model = setParam(model, 'lb', objectiveFunction, atpRates(i));
        model = setParam(model, 'ub', objectiveFunction, atpRates(i)*1.001);
        
        if settings.pfba == true
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
