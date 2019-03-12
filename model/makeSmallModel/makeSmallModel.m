%load('reducedModel')
load('../muscleModel')
addpath('../../src1')

objectiveFunction = {'human_ATPMaintainance'};
ATPProduction = findIndex(model.rxns, objectiveFunction);

%Minimal media simulation
inMedia = {
    'octanoic acid[s]'
    'glycogen[s]'
    'O2[s]'  
    'H2O[s]'
};

outMedia = {
    'CO2[s]' 
    'H2O[s]'
    'L-lactate[s]'
};



[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

%Remove ProtPool

reactionIn = getBounds(model, inMedia);
reactionOut = getBounds(model, outMedia);

model = setParam(model, 'lb', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', exchangeRxnsIndexes, 0);
model = setParam(model, 'ub', reactionOut, 1000);
model = setParam(model, 'lb', reactionIn, -1000);
model = setParam(model, 'lb', reactionIn(1), 0);
model = setParam(model, 'lb', reactionIn(2), 0);

model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);

% banSubsystems = {
%     'Bile acid biosynthesis'
%     'Glutathione metabolism'
%     'Porphyrin metabolism'
%     'Eicosanoid metabolism'
%     'Serotonin and melatonin biosynthesis'
%     'Folate metabolism'
%     };
% model = removeSubsystems(model, banSubsystems);


exceptions = {'HMR_3949', 'HMR_3825','HMR_4852','HMR_4888','HMR_4898','HMR_4922','HMR_4926','HMR_5043','HMR_6328','HMR_7638'};
mitTransport = ismember(model.subSystems, 'Transport, mitochondrial');
for i = 1:length(exceptions)
    mitTransport(findIndex(model.rxns, exceptions{i})) = 0;
end
model.lb(mitTransport) = 0;
model.ub(mitTransport) = 0;

%Rename subsystems:
model.subSystems{findIndex(model.rxns, 'HMR_4137')} = 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism';
model.subSystems(ismember(model.subSystems, 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism')) = {'Tricarboxylic acid cycle'};



%Remove all compartments appart from cytosol, mitochondria and extracell
model.lb(ismember(model.rxnComps, [2 5 6 7 8])) = 0;
model.ub(ismember(model.rxnComps, [2 5 6 7 8])) = 0;

%Prevent cytsolic TCA
% model.lb(findIndex(model.rxns,'HMR_4957')) = 0;
% model.ub(findIndex(model.rxns,'HMR_4957')) = 0;




%Fatty acid flux
model = setParam(model, 'lb', reactionIn(1), -1);
solution0 = solveLinMin(model, true)


%Canonical flux
model = setParam(model, 'lb', reactionIn(1), 0);
model = setParam(model, 'lb', reactionIn(2), -1);
solutionA = solveLinMin(model, true)

printAffected(model, solutionA, solution0)

%Complex 1 knock out
model1 = setParam(model, 'lb', 'HMR_6921', 0);
model1 = setParam(model1, 'ub', 'HMR_6921', 0);

solutionB = solveLinMin(model1, true)
%printAffected(model, solutionA, solutionB);


%glycerol-3-phosphate shunt knock
model2 = setParam(model1, 'lb', 'HMR_0483', 0);
model2 = setParam(model2, 'ub', 'HMR_0483', 0);
solutionC = solveLinMin(model2, true)
%printAffected(model, solutionB, solutionC)




%ETF-ubiquinone oxidoreductase knock
% model = setParam(model, 'lb', 'HMR_6911', 0);
% model = setParam(model, 'ub', 'HMR_6911', 0);

% solutionD = solveLinMin(model, true)
% printAffected(model, solutionC, solutionD)

%ATP synthase KO
model3 = setParam(model, 'lb', 'HMR_6916', 0);
model3 = setParam(model3, 'ub', 'HMR_6916', 0);

solutionD = solveLinMin(model3, true)
%printAffected(model, solutionC, solutionD)

%oxygen transport
model4 = setParam(model, 'lb', 'HMR_9048', 0);
model4 = setParam(model4, 'ub', 'HMR_9048', 0);
solutionE = solveLinMin(model4, true)

totalSolution = abs(solution0.x) + abs(solutionA.x) + abs(solutionB.x) + abs(solutionC.x) + abs(solutionD.x) + abs(solutionE.x);

% clc
% printExchangeFluxes(model, solutionE.x)
tresh = 10^-5;
reducedModel = removeRxns(model,totalSolution<tresh,true,true,true);

model = reducedModel;
save('reducedModel.mat', 'model');

reducedModel = rmfield(reducedModel,'rxnComps');
reducedModel = rmfield(reducedModel,'geneComps');
exportToExcelFormat(reducedModel,'../reducedModel.xls')
