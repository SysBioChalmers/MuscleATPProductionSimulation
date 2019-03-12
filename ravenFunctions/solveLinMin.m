function solution = solveLinMin(model, silent, epsilon)

    if nargin <3
        epsilon = 0;
    end
    
    if nargin <2
       silent = false; 
    end
        
    nrOfReactions = length(model.c);   
    
    firstSolution = solveLin(model, silent);
    
    if length(firstSolution.x) == nrOfReactions
        objective = find(model.c);
        objectiveValues = firstSolution.x(objective);
        
        objectivePositive = objective(objectiveValues>=0);
        objectivePositiveValues = objectiveValues(objectiveValues>=0);
        objectiveNegative = objective(objectiveValues<=0);
        objectiveNegativeValues = -objectiveValues(objectiveValues<=0);

        %Add reversed Matrix
        newModel = model;
        newModel.S = horzcat(model.S, -model.S);    
        newModel.lb = vertcat(model.lb,  -model.ub);
        newModel.ub = vertcat(model.ub, -model.lb);
        newModel.c = vertcat(model.c, -model.c);

        %Add objectiveFlux
        newModel.lb(objectivePositive) = objectivePositiveValues - epsilon;
        newModel.ub(objectivePositive) = objectivePositiveValues + epsilon;
        newModel.ub(objectiveNegative) = 0;
        newModel.lb(nrOfReactions + objectiveNegative) = objectiveNegativeValues + epsilon;
        newModel.ub(nrOfReactions + objectiveNegative) = objectiveNegativeValues - epsilon;
        
        newModel.ub(nrOfReactions + objectivePositive) = 0;


        %Set Lower Bound
        newModel.lb(newModel.lb<0) = 0;
        newModel.ub(newModel.ub<0) = 0;

        %Minimize fluxes
        newModel.c = -ones(nrOfReactions*2,1);

        solution = solveLin(newModel, true);
        
        if epsilon == 0 && length(solution.x) < nrOfReactions
            disp('Warning, no solution found perturbing arround solution with an epsilon of 10^-5')
            testEpsilon = 10^-5;
            solution = solveLinMin(model, silent, testEpsilon);
        elseif length(solution.x) < nrOfReactions
            disp('Warning, could not solve second optimization returning first')
            solution = firstSolution;           
        else
            solution.x = solution.x(1:nrOfReactions) - solution.x((nrOfReactions+1):end);
            solution.f = sum(solution.x(objectivePositive)) - sum(solution.x(objectiveNegative));
            solution.f = -solution.f;
        end
    else
        solution = firstSolution;
    end
end