load('../muscleModel');

%Add ATP work "metabolite"
glycogenPhos = createRXNStuct(model, 'ATPwork', 'ATP[c] + H2O[c] => ADP[c] + Pi[c] + ATPwork[c]', 0, 1000, 'Artificial');
model=addRxns(model,glycogenPhos,3,'c',true);
glycogenPhos = createRXNStuct(model, 'ATPworkOut', 'ATPwork[c] => ATPwork[s]', 0, 1000, 'Artificial');
model=addRxns(model,glycogenPhos,3,'c',true);

%Muscle Model
modelName = 'MuscleConnect';
ids = {'m1', 'm2', 'm3'};
names = {'Muscle1', 'Muscle2', 'Muscle3'};
models{1} = model;
models{2} = model;
models{3} = model;

compartmentName = 'blood';
compartmentId = 'b';
superModel = buildModel(models, modelName, ids, names, compartmentId, compartmentName);

superModel = addCompartment(superModel, '', 's');


[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(superModel,'both');
constructEquations(superModel, exchangeRxnsIndexes)


superModel.compNames = superModel.comps;

superModel = addSupportReactions(superModel, 0.5);

%Add lactate accumulation
%(Lactate- H+) + (HCO3- K+) -> (K+ Lactate-) + CO2 + H2O
lactateBuffer = createRXNStuct(superModel, 'lactateBuffering', 'L-lactate[sb] => CO2[sb] + H2O[sb]', 0, 1000, 'Artificial reactions');
superModel=addRxns(superModel,lactateBuffer,3,'sb',false);

save('connectedMuscles.mat', 'superModel')