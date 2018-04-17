function [ output_args ] = plotRelevantPlots(model, x, y, exchangeRxnsIndexes, fat, oxygen, co2)
    rxns = constructEquations(model, exchangeRxnsIndexes);
    %remove crap
    for i = 1:length(rxns)
        tmp = strsplit(rxns{i}, ' <=>');
        rxns{i} = tmp{1}; 
    end
    fatRxn = ismember(rxns, fat);
    oxygenRxn = ismember(rxns, oxygen);
    co2Rxn = ismember(rxns, co2);
    
    vo2 = abs(y(:,oxygenRxn));
    vo2Max = vo2./max(vo2);
    co2 =  abs(y(:,co2Rxn));
    RER = co2./vo2;
    
    plot(vo2Max, -y(:,fatRxn));
    figure()
    plot(vo2Max, RER);
    
end

