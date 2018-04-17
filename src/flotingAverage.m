function [growthRates, output ] = flotingAverage(model, growthRates, results, range)
    compForFlux = findIndex(model.comps, 'sb');
    biomassName = 'human_biomass';
    fluxThresh = 10^-3;
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
    
   
    hold all       
    values = zeros(size(model.S, 1), size(results,1));
    for j = 1:size(results,1)
        values(:,j) = model.S(:,exchangeRxnsIndexes) * results(j, exchangeRxnsIndexes)';
    end

    metabolites = abs(sum(values, 2))>fluxThresh;
    metNames = model.metNames(metabolites);
    metComp = model.metComps(metabolites);
    metNames = metNames(metComp == compForFlux);
    values = values(metabolites,:);
    values = values(metComp == compForFlux,:);
    biomassMet = ismember(metNames, biomassName);
    metNames(biomassMet) = [];
    values(biomassMet,:) = [];

    metNames = indicateDirection(metNames, values);
    
    lag = round(length(growthRates)*range);

    testGrowth = tsmovavg(growthRates,'s',lag,2);
    
    output = tsmovavg(values','s',lag,1);
    plot(testGrowth, abs(output), 'linewidth', 3)
    legend(metNames, 'location', 'nw');
end

function metNames = indicateDirection(metNames, values)
    directionality = sign(mean(values,2));
    for i =1:length(metNames)
        if directionality(i) == 1
            metNames{i} = ['=>' metNames{i}];
        else
            metNames{i} = ['<=' metNames{i}];
        end
    end
end
