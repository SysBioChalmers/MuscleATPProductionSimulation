function reaction = getTransportReactions(model, comp1, comp2)
    compartmentNr1 = findIndex(model.comps, comp1);
    compartmentNr2 = findIndex(model.comps, comp2);
    compList1 = ismember(model.metComps, compartmentNr1);
    compList2 = ismember(model.metComps, compartmentNr2);
        
    reaction1 = sum(abs(model.S(compList1,:)),1)>0;
    reaction2 = sum(abs(model.S(compList2,:)),1)>0;
    
    reaction = and(reaction1, reaction2);
end

