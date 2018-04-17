function [exchangeRxnsIndexes, results] = plotFullSolution(model, growthRates, fullSolution, objectiveFunction)
    fluxThresh = 10^-3;
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
    
    skipBiomass = findIndex(model.rxns, objectiveFunction);
    
    exchangeRxnsIndexes(exchangeRxnsIndexes==skipBiomass) = [];  
    
    results = fullSolution(:,exchangeRxnsIndexes);

    
    hold all

    highestFlux = max(max(abs(results)));


    emptySolutions = mean(abs(results),2)<=0;
    results(emptySolutions,:) = [];
    growthRates(emptySolutions) = [];
   
    
    if not(isempty(growthRates))

        emptyExchange = sum(abs(results))<= (fluxThresh * highestFlux);

        results(:, emptyExchange) = [];
        exchangeRxnsIndexes(emptyExchange) = [];

        positiveFluxes = mean(results)>0;

        
        for i=0:1
            subplot(2,1,i+1);
            plot(growthRates, abs(results(:,positiveFluxes == i)) , 'linewidth', 3)
            eqns=constructEquations(model, model.rxns(exchangeRxnsIndexes(positiveFluxes == i)));
            legend(eqns, 'location', 'nw')
            ylim([0 highestFlux]);
            xlim([min(growthRates) max(growthRates)]);
            xlabel('Rate', 'FontSize',15,'FontName', 'Arial')
            ylabel('Flux mMol/gdw/h', 'FontSize',15,'FontName', 'Arial')
            set(gca,'FontSize',15,'FontName', 'Arial')
        end

    hold off


end