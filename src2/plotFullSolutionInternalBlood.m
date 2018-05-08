function data = plotFullSolutionInternalBlood(model, growthRates, results, metaboliteList, bloodComp, otherComp)
    hold all       
    targetValue = 0.3;
    ylimValue = 7;
    ytext = 'Flux mol/h';    
    xtext = 'Watt';
    W = 1000*molToW(growthRates);
    
        for i=1:length(otherComp)
            subplot(1,length(otherComp),i);
            
            

%           [metabolites, values] = getAllMetabolites(model, results, reactionSet{i});
            transportReactions = getTransport(model, metaboliteList, bloodComp, otherComp{i});
            values = results(:,transportReactions)';

            metNames = metaboliteList;           
            metNames = indicateDirection(metNames, values);
            
%            targetValue = max(max(values)) * amplificationThresh;
            
            for j = 1:size(values,1)
               if max(abs(values(j,:)))< targetValue 
                    values(j,:) = values(j,:)*10;
                    metNames{j} = [metNames{j} '[x10]'];
               end
            end
            
            %Show legend and y axis in subfigure 1
            if i>1
                ytext = [];
            end

            if i<length(otherComp)
                metNames = [];
            end
            
            
            nicePlot(W, abs(values'), metNames, xtext, ytext, otherComp{i});
            ylim([0 ylimValue]);
        end
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

function [metabolites, values] = getAllMetabolites(model, results, currentReactions)
    fluxThresh = 10^-3;
    values = zeros(size(model.S, 1), size(results,1));
    for j = 1:size(results,1)
        values(:,j) = model.S(:,currentReactions) * results(j, currentReactions)';
    end
    metabolites = abs(sum(values, 2))>fluxThresh;
end


function nicePlot(x, y, ltext, xtext, ytext, ttext)
    plot(x, y, 'linewidth', 2)

    if not(isempty(ltext))
        legend(ltext, 'location', 'nw')
        legend boxoff
    end
    xlim([min(x) max(x)]);
    xlabel(xtext, 'FontSize',12,'FontName', 'Arial')
    if not(isempty(ytext))
        ylabel(ytext, 'FontSize',12,'FontName', 'Arial')
    end
    set(gca,'FontSize',12,'FontName', 'Arial')
    title(ttext)
end