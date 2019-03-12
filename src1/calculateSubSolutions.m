function [glyYields, massYield, solutionNames, solutions] = calculateSubSolutions(model)
    massConstraintRow = findIndex(model.mets, 'MassConstraint');
    model.b(massConstraintRow,2) = 1000; %unconstrained mass for these solutions
    weightRow = full(model.S(massConstraintRow,:));

    minimalMedia = {
        'O2[s]'
        'octanoic acid[s]'
        'glycogen[s]'
        'H2O[s]'
        };

    minimalFlux= [-1000
                  0
                  0
                  -1000];

    model = configureModel(model, minimalMedia, minimalFlux);
    
    reactionNumbers = getBounds(model, minimalMedia);
    fatNr = reactionNumbers(2);
    glucoseNr = reactionNumbers(3);

    objectiveFunction = {'human_ATPMaintainance'};

    model = setParam(model, 'lb', objectiveFunction, 0);
    model = setParam(model, 'ub', objectiveFunction, 1000);
    model = setParam(model, 'obj', objectiveFunction, 1);

    %solutionNames = {'Aerobic', 'Complex I bypass', 'Uncoupling', 'Fermentative'};
    solutionNames = {'Aerobic (fat)', 'Aerobic (glycogen)', 'Complex I bypass', 'Fermentative', 'Uncoupling'};
    solutions = zeros(5,length(model.rxns));

    %Fatty acid:
    model = setParam(model, 'lb', fatNr, -1); 
    res = solveLinMin(model);
    solutions(1,:) = res.x;
    
    %Remove lactate
    model = setParam(model, 'lb', fatNr, 0); 
    model = setParam(model, 'lb', glucoseNr, -1); 
    res = solveLinMin(model);
    solutions(2,:) = res.x;


    %bypass of complex 1
    model = setParam(model, 'ub', 'HMR_6921', 0); 
    res = solveLinMin(model);
    solutions(3,:) = res.x;

    %uncoupling
    model = setParam(model, 'ub', 'HMR_6916', 0); 
    res = solveLinMin(model);
    solutions(5,:) = res.x;

    %Remove oxygen
    model = setParam(model, 'ub', 'HMR_6914', 0); 
    res = solveLinMin(model);
    solutions(4,:) = res.x;

    fullMass = abs(solutions) .* repmat(weightRow, size(solutions,1),1);
    massOfFlux = sum(fullMass,2);

    cMol = -6*solutions(:,glucoseNr)' + -8*solutions(:,fatNr)';
    
    atpFlux = solutions(:,findIndex(model.rxns, objectiveFunction));

    glyYields = atpFlux./cMol';
    massYield = atpFlux./massOfFlux;


    solutions = solutions./repmat(atpFlux, 1, size(solutions,2));
end

