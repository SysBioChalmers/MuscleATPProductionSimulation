function printAffected(model, solutionA, solutionB)

    tresh = 10^-5;
    rxnsWithFluxA = find(abs(solutionA.x)>tresh);
    rxnsWithFluxB = find(abs(solutionB.x)>tresh);
    
    [C,ia,ib] = setxor(rxnsWithFluxA,rxnsWithFluxB);
    eqA = constructEquations(model, rxnsWithFluxA(ia));
    eqB = constructEquations(model, rxnsWithFluxB(ib));
    
    rxnA = model.rxns(rxnsWithFluxA(ia));
    rxnB = model.rxns(rxnsWithFluxB(ib));
    
    for i = 1:length(eqA)
        fprintf('%s\t%s\n', rxnA{i}, eqA{i})
    end
    fprintf('\n\n')
    
    for i = 1:length(eqB)
        fprintf('%s\t%s\n', rxnB{i}, eqB{i})
    end    
    
end

