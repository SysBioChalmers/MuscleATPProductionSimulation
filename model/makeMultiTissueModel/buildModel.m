function superModel =  buildModel(subModels, modelName, subId, subNames, compId, compartmentName)
    superModel.id = modelName;
    
    superModel.subModels = subModels;
    for i = 1:length(subModels)
        superModel.subModels{i}.id = subId{i};
        superModel.subModels{i}.name = subNames{i};
    end   

    superModel = mergeModels(superModel);    
    superModel = addCompartment(superModel, compId, compartmentName);
end

