function superModel = mergeModels(superModel)
    curPos = [0 0];
    superModel = initStruct(superModel);
    
    allSize = zeros(length(superModel.subModels),2);
    for i = 1:length(superModel.subModels)
        allSize(i, :) = size(superModel.subModels{i}.S);
    end
    superModel.S = sparse(sum(allSize));
    
    for i = 1:length(superModel.subModels)
        curModel = superModel.subModels{i};
        rInterval = (1 + curPos(1)):(curPos(1)+allSize(i, 1));
        cInterval = (1 + curPos(2)):(curPos(2)+allSize(i, 2));
        superModel.S(rInterval, cInterval) = curModel.S;
        curPos = curPos + allSize(i,:);
        
        superModel = appendStruct(superModel, curModel);
    end 
end




function superModel = initStruct(superModel)
    superModel.rxns = [];
    superModel.mets = [];
    superModel.rxnNames = [];
    superModel.metNames = [];
    %superModel.rxnComps = [];
    superModel.subSystems = [];
    superModel.metComps = [];  
    superModel.comps = [];
    superModel.lb = [];
    superModel.ub = [];
    superModel.rev = [];
    superModel.b = [];
    superModel.c = [];
end
