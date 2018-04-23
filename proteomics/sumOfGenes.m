function subMass = sumOfGenes(geneMap, involvedGenes)
    subMass = 0;
    for j = 1:length(involvedGenes)
       if isKey(geneMap, involvedGenes{j})
          subMass = subMass + geneMap(involvedGenes{j});
       else
          disp(involvedGenes{j})
       end
    end
end

