function [ output_args ] = printFluxesValues(model, results)
    nonZero = sum(abs(results))>10^-6;
    fluxes = results(:, nonZero);
    subSystem = model.subSystems(nonZero);
    rxns = model.rxns(nonZero);
    eq = constructEquations(model, model.rxns(nonZero));
    
    
    [values, order] = sort(sum(abs(fluxes)));
    
    for i = 1:length(values)
        curRow = order(i);
        organ = strsplit(rxns{curRow}, '_');
        organ = organ{end};
        fprintf('%s\t%s\t%s\t%s', rxns{curRow}, subSystem{curRow}, organ, eq{curRow}); 
        for j=1:size(fluxes,1)
            fprintf('\t%f', fluxes(j, curRow)); 
        end
        fprintf('\n'); 
    end
    
    
end

