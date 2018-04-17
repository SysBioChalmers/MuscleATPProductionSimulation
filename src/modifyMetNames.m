function metNames = modifyMetNames(modelMod)
    metNames = modelMod.metNames;
    for i=1:length(metNames)
        metNames{i} = [metNames{i} '[' modelMod.compNames{modelMod.metComps(i)} ']'];
    end
end

