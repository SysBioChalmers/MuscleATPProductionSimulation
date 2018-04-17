function data = plotFullSolutionAandB(model, growthRates, results, bloodComp, AandB)
    compForFlux = findIndex(model.comps, bloodComp);
    skipList = {'obectiveMetabolite', 'bloodMet', 'H2O'};
    fluxThresh = 10^-3;

    amplificationThresh = 0.05;

    reactionSet{1} = find(getTransportReactions(model, bloodComp, AandB{1}))';
    reactionSet{2} = find(getTransportReactions(model, bloodComp, AandB{2}))';

    hold all       
    valuesA = zeros(size(model.S, 1), size(results,1));
    valuesB = zeros(size(model.S, 1), size(results,1));
    for j = 1:size(results,1)
        valuesA(:,j) = model.S(:,reactionSet{1}) * results(j, reactionSet{1})';
        valuesB(:,j) = model.S(:,reactionSet{2}) * results(j, reactionSet{2})';
    end
    values = valuesA + valuesB;
        
    metabolites = abs(sum(values, 2))>fluxThresh;
    metIndex = find(metabolites);
    metNames = model.metNames(metIndex);
    metComp = model.metComps(metIndex);
    metIndex = metIndex(metComp == compForFlux);
    metNames = metNames(metComp == compForFlux);
    metIndex = metIndex(ismember(metNames, skipList)==0);
    metNames = metNames(ismember(metNames, skipList)==0);
       
    values = values(metIndex,:);
    lactrxn = find(ismember(metNames,'L-lactate'));
    O2rxn = find(ismember(metNames,'O2'));
    data = [1000*molToW(growthRates)', -molToMl(values(O2rxn,:))', 1000/60*values(lactrxn,:)'];
    %    for i = 1:length(data)
    %       fprintf('%f\t%f\t%f\n', data(i,1), data(i,2), data(i,3));
    %    end

    
    metNames = indicateDirection(metNames, values);

    targetValue = max(max(values)) * amplificationThresh;

    for j = 1:size(values,1)
       if max(abs(values(j,:)))< targetValue 
            values(j,:) = values(j,:)*10;
            metNames{j} = [metNames{j} '[x10]'];
       end
    end


    nicePlot(growthRates, abs(values'), metNames, 'Rate', 'Flux mol/h', 'A+B');

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