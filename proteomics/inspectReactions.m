function inspectReactions(model, allRxns, geneMap, tresh)
results = zeros(length(allRxns),1);
allEqns = constructEquations(model, allRxns);
for i = 1:length(allRxns)
   curRxns = ismember(model.rxns, allRxns{i});
   involvedGenes = getInvolvedGenes(model, curRxns);
   subMass = 0;
   for j = 1:length(involvedGenes)
       if isKey(geneMap, involvedGenes{j})
          subMass = subMass + geneMap(involvedGenes{j});
       else
          disp(involvedGenes{j})
       end
   end
   results(i) = subMass;
end

%normalize
results = results/sum(results);

figure()
[results, indx] = sort(results, 'descend');
allEqns = allEqns(indx);
allEqns(results<tresh) = [];
results(results<tresh) = [];
figure()
plotPie(results, allEqns)
allEqns
end

