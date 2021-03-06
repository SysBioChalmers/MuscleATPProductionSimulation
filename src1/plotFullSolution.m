function plotFullSolution(model, growthRates, fullSolution, plotExchange)
    fluxThresh = 10^-3;
 
    
    reactionNumbers= getBounds(model, plotExchange);
    plotExchange = simplifyMetNames(plotExchange);
    
    
    results = fullSolution(:,reactionNumbers);
 
   
    hold all

    highestFlux = max(max(abs(results)));


%     emptySolutions = mean(abs(results),2)<=0;
%     results(emptySolutions,:) = [];
%     growthRates(emptySolutions) = [];
    
    
    if not(isempty(growthRates))
        
        maxYVal = ceil(max(abs(results(:)))*1.1);
        maxXval = max(growthRates)*1.1;

        positiveFluxes = mean(results)>0;
        results(:,not(positiveFluxes)) =  -results(:,not(positiveFluxes));
        
        plot(growthRates, results, 'linewidth', 3)

        legend(plotExchange, 'location', 'nw')
        legend boxoff
        ylim([0 maxYVal]);
        xlim([0 maxXval]);
        xlabel('Rate', 'FontSize',14,'FontName', 'Arial')
        ylabel('Flux mmol/gdw/h', 'FontSize',14,'FontName', 'Arial')
        set(gca,'FontSize',14,'FontName', 'Arial')
        %figure()
        %plot(growthRates, weightEstimate)
    end



end

function eqn2 = getMetabolites(model, rxns)
    eqn2=constructEquations(model, model.rxns(rxns));
    for i = 1:length(eqn2)
       eqn2{i} = strrep(eqn2{i}, '[s]', '');
       eqn2{i} = strrep(eqn2{i}, '<=>', '');       
    end
end

function mets = simplifyMetNames(mets)
    for i = 1:length(mets)
       mets{i} = strrep(mets{i}, '[s]', '');      
    end
end