function plotFullSolution(model, growthRates, fullSolution, objectiveFunction)
    fluxThresh = 10^-3;
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
    
    skipBiomass = findIndex(model.rxns, objectiveFunction);
    
    exchangeRxnsIndexes(exchangeRxnsIndexes==skipBiomass) = []; 
    
    %remove water
    exchangeRxnsIndexes(exchangeRxnsIndexes==getBounds(model, {'H2O[s]'})) = []; 
    
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
        
        maxYVal = max(abs(results(:)))*1.1;
        maxXval = max(growthRates)*1.1;

        positiveFluxes = mean(results)>0;
        
        eqn1=getMetabolites(model, exchangeRxnsIndexes(not(positiveFluxes)));
        eqn2=getMetabolites(model, exchangeRxnsIndexes(positiveFluxes));
        
        subplot(2,1,1);
        plot(growthRates, -results(:,not(positiveFluxes)) , 'linewidth', 3)

        legend(eqn1, 'location', 'nw')
        legend boxoff
        set(gca,'FontSize',14,'FontName', 'Arial')
        set(gca,'XTick',[])
        ylim([0 maxYVal]);
        xlim([0 maxXval]);
        subplot(2,1,2);
        plot(growthRates, results(:,positiveFluxes), 'linewidth', 3)

        
        model.rxns(exchangeRxnsIndexes(positiveFluxes))
        legend(eqn2, 'location', 'nw')
        legend boxoff
        ylim([0 maxYVal]);
        xlim([0 maxXval]);
        xlabel('Rate', 'FontSize',14,'FontName', 'Arial')
        ylabel('Flux mMol/gdw/h', 'FontSize',14,'FontName', 'Arial')
        set(gca,'FontSize',14,'FontName', 'Arial')
        %figure()
        %plot(growthRates, weightEstimate)
    end

    hold off


end

function eqn2 = getMetabolites(model, rxns)
    eqn2=constructEquations(model, model.rxns(rxns));
    for i = 1:length(eqn2)
       eqn2{i} = strrep(eqn2{i}, '[s]', '');
       eqn2{i} = strrep(eqn2{i}, '<=>', '');       
    end
end