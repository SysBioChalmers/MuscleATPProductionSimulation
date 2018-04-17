function [exchange, values] = plotFullSolutionInternal(model, growthRates, results, objectiveFunction)
    fluxThresh = 10^-3;
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');
    transportReactions1 = find(getTransportReactions(model, 'sm', 'sb'));
    transportReactions2 = find(getTransportReactions(model, 'sb', 'sl'));
    
    
    skipBiomass = findIndex(model.rxns, objectiveFunction);
    exchangeRxnsIndexes(exchangeRxnsIndexes==skipBiomass) = [];  
    reactionSet{1} = exchangeRxnsIndexes;
    reactionSet{2} = transportReactions1;
    reactionSet{3} = transportReactions2;
    
    titles = {
        'Exchange'
        'B->M'
        'B->L'
        };
    
    hold all       
        for i=1:3
            subplot(2,2,i);
            currentReactions = reactionSet{i};
            emptyReactions = mean(abs(results(:,currentReactions)))<10^-5;
            currentReactions(emptyReactions) = [];
            plot(growthRates, abs(results(:,currentReactions)), 'linewidth', 3)
            eqns = cleanEquations(model, results, currentReactions);
            legend(eqns, 'location', 'nw')
            %ylim([0 highestFlux]);
            xlim([min(growthRates) max(growthRates)]);
            xlabel('Rate', 'FontSize',15,'FontName', 'Arial')
            ylabel('Flux mMol/gdw/h', 'FontSize',15,'FontName', 'Arial')
            set(gca,'FontSize',15,'FontName', 'Arial')
            title(titles{i})
            
            if i ==1
                exchange = currentReactions;
                values = results(:,currentReactions);
            end
        end

    hold off
end


function eqns = cleanEquations(model, results, currentReactions)
    data = results(:, currentReactions);
    directionality = sign(mean(data));
    eqns = constructEquations(model, model.rxns(currentReactions));
    for i =1:length(eqns)
        tmp = strsplit(eqns{i}, ' <=>');
        eqns{i} =  tmp{1};
        if directionality(i) == 1
            eqns{i} = ['=>' eqns{i}];
        else
            eqns{i} = ['<=' eqns{i}];
        end
    end

end