function indx = getComplexIbypass(model, fullSolution, rxn)
    tresh = 10^-6;
    complexIvalues = fullSolution(:,findIndex(model.rxns, rxn));
    complexIvalues = complexIvalues/max(complexIvalues);
    indx = find(complexIvalues > 1-tresh);
    indx = indx(1);
end

