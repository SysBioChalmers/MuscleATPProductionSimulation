function affectedReactions = getReactionsFromCompartment(model, comp)
    compartmentNr = findIndex(model.comps, comp);
    affectedMets = ismember(model.metComps, compartmentNr);
    affectedReactions = sum(abs(model.S(affectedMets, :)),1)>0;
end

