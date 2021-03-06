function [glyYields, massYield, solutionNames, solutions] = calculateSubSolutions(model)
    massConstraintRow = findIndex(model.mets, 'MassConstraint');
    model.b(massConstraintRow,2) = 1000; %unconstrained mass for these calculations
    weightRow = full(model.S(massConstraintRow,:));

    minimalMedia = {
        'O2[s]'
        'palmitate[s]'
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
    
    model = setParam(model, 'lb', fatNr, 0);
    model = setParam(model, 'lb', glucoseNr, 0);
    

    objectiveFunction = {'human_ATPMaintainance'};

    model = setParam(model, 'lb', objectiveFunction, 0);
    model = setParam(model, 'ub', objectiveFunction, 1000);
    model = setParam(model, 'obj', objectiveFunction, 1);

    solutionNames = {'oxidative (fat)', 'oxidative (glycosyl)',  'oxidative (glycosyl)', 'complex I bypass', 'fermentative', 'uncoupling'};
    solutions = zeros(length(solutionNames),length(model.rxns));

    tresh = 10^-6;
    for i = 1:size(solutions,1)
        tmpModel = model;
        switch i
            case 1
                %Fatty acid:
                tmpModel = setParam(tmpModel, 'lb', fatNr, -1);
            case 2
                %Glucose
                tmpModel = setParam(tmpModel, 'lb', glucoseNr, -1);
            case 3
                %Glucose
                tmpModel = setParam(tmpModel, 'lb', glucoseNr, -1);
                %Block Mal-Asp
                tmpModel = setParam(tmpModel, 'ub', 'HMR_3829', 0);
            case 4
                %bypass of complex 1
                tmpModel = setParam(tmpModel, 'lb', glucoseNr, -1);
                tmpModel = setParam(tmpModel, 'ub', 'HMR_6921', 0);
            case 5              
                %Remove oxygen
                tmpModel = setParam(tmpModel, 'lb', glucoseNr, -1);
                tmpModel = setParam(tmpModel, 'ub', 'HMR_6914', 0);               
            case 6
                %uncoupling
                tmpModel = setParam(tmpModel, 'lb', glucoseNr, -1);
                tmpModel = setParam(tmpModel, 'ub', 'HMR_6916', 0);                 
        end
        
        %Minimize protein instead of flux
        %res = solveLinMin(tmpModel);
        res = solveLin(tmpModel);
        tmpModel.lb(find(tmpModel.c)) = -res.f -tresh;
        tmpModel.c = -weightRow;
        res = solveLin(tmpModel);
        
        solutions(i,:) = res.x;
        %printFluxes(tmpModel, res.x, false, 0, [], '%rxnID\t%eqn\t%flux\n')
    end

    fullMass = abs(solutions) .* repmat(weightRow, size(solutions,1),1);
    massOfFlux = sum(fullMass,2);

    cMol = -6*solutions(:,glucoseNr)' + -16*solutions(:,fatNr)';
    
    atpFlux = solutions(:,findIndex(model.rxns, objectiveFunction));

    glyYields = atpFlux./cMol';
    massYield = atpFlux./massOfFlux;

    solutions = solutions./repmat(atpFlux, 1, size(solutions,2));
end

