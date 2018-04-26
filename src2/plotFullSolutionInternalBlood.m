function data = plotFullSolutionInternalBlood(model, growthRates, results, bloodComp, otherComp)

    compForFlux = findIndex(model.comps, bloodComp);
    skipList = {'obectiveMetabolite', 'bloodMet', 'H2O', 'bloodM1', 'bloodM2'};
    fluxThresh = 10^-3;
    [exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

    amplificationThresh = 0.05;

    reactionSet{1} = exchangeRxnsIndexes;
    reactionSet{1} = find(getTransportReactions(model, bloodComp, 's'));
    titles{1} = 'Exchange';
    for i = 1:length(otherComp)
        transportReactions = find(getTransportReactions(model, bloodComp, otherComp{i}))';
        reactionSet{i+1} = transportReactions;
        titles{i+1} = otherComp{i};
    end
        
    hold all       
        for i=1:length(reactionSet)
            subplot(2,2,i);
            currentReactions = reactionSet{i};
            values = zeros(size(model.S, 1), size(results,1));
            for j = 1:size(results,1)
                values(:,j) = model.S(:,currentReactions) * results(j, currentReactions)';
            end
            
            metabolites = abs(sum(values, 2))>fluxThresh;
            metNames = model.metNames(metabolites);
            metComp = model.metComps(metabolites);
            metNames = metNames(metComp == compForFlux);
            values = values(metabolites,:);
            values = values(metComp == compForFlux,:);
            biomassMet = ismember(metNames, skipList);
            metNames(biomassMet) = [];
            values(biomassMet,:) = [];
            
            metNames = indicateDirection(metNames, values);
            
            targetValue = max(max(values)) * amplificationThresh;
            
            for j = 1:size(values,1)
               if max(abs(values(j,:)))< targetValue 
                    values(j,:) = values(j,:)*10;
                    metNames{j} = [metNames{j} '[x10]'];
               end
            end
            %growthRates'
            %abs(values')
            data{i} = abs(values');
            nicePlot(growthRates, abs(values'), metNames, 'Rate', 'Flux mol/h', titles{i});
        end
%     subplot(2,2,length(reactionSet)+1);
%     
%     
%     bloodReaction= findIndex(model.rxns, 'BloodCost');
%     if ~isempty(bloodReaction)
%         bloodValues = results(:,bloodReaction);
%         bloodMets = find(model.S(:,bloodReaction));
%         bioBlood = ismember(model.metNames(bloodMets), skipList);
%         bioBlood = bloodMets(bioBlood);
%         bloodFactor = model.S(bioBlood,bloodReaction);
%         if isempty(bloodFactor)
%             bloodFactor = 0;
%         end        
%         
%         ratio = -bloodFactor(1) * full(bloodValues)./growthRates';
%         linFit = polyfit(growthRates(1:(end/2)), bloodValues(1:(end/2))',1);
%         linReg = linFit(2) + linFit(1).*growthRates;
%         legendValues = {'Blood Flow', 'ATPCost', 'linReg BlodFLow'};
%         nicePlot(growthRates, [bloodValues ratio linReg'], legendValues, 'Rate', 'Unit', 'BloodFlow');    
%     
%         lactateMet = ismember(model.metNames, 'L-lactate');
%         liverTransport = find(getTransportReactions(model, bloodComp, 'sl'));
%         lactateRxn = find(sum(abs(model.S(lactateMet, liverTransport)),1));
%         lactateRxn = liverTransport(lactateRxn);
%         lactateFlux = -results(:,lactateRxn);
% 
%         CO2Met = ismember(model.metNames, 'CO2');
%         CO2Rxn = find(sum(abs(model.S(CO2Met, exchangeRxnsIndexes)),1));
%         CO2Rxn = exchangeRxnsIndexes(CO2Rxn);
%         cO2Flux = results(:,CO2Rxn);    
% 
%         subplot(2,2,length(reactionSet)+2);
%         legendValues = {'lactate', 'CO2'};
%         lactateSaturation = lactateFlux./full(bloodValues);
%         cO2FluxSaturation = cO2Flux./full(bloodValues);
%         nicePlot(growthRates, [lactateSaturation cO2FluxSaturation], legendValues, 'Rate', 'Unit', 'Saturation');         
%         
%         
%     end
 
    
    hold off
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


function nicePlot(x, y, ltext, xtext, ytext, ttext)
    plot(x, y, 'linewidth', 3)

    if not(isempty(ltext))
        legend(ltext, 'location', 'nw')
    end
    xlim([min(x) max(x)]);
    xlabel(xtext, 'FontSize',12,'FontName', 'Arial')
    ylabel(ytext, 'FontSize',12,'FontName', 'Arial')
    set(gca,'FontSize',12,'FontName', 'Arial')
    title(ttext)
end