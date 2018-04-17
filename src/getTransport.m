function result = getTransport(model, metabolite, comp1, comp2)
    allReactions = getTransportReactions(model, comp1, comp2);

    result = zeros(length(metabolite),1);
    reducedS = abs(model.S(:,allReactions));
    allReactonNr = find(allReactions);
    
    for i = 1:length(metabolite)
       currentMets = ismember(model.metNames, metabolite{i});
       currentReaction = find(sum(reducedS(currentMets,:),1));
       allReactonNr(currentReaction);
       result(i) = allReactonNr(currentReaction);
    end
end

