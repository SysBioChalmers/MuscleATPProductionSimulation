%load('reducedModel')
load('../muscleModel')
addpath('../../src1')

objectiveFunction = {'human_ATPMaintainance'};
ATPProduction = findIndex(model.rxns, objectiveFunction);

%Minimal media simulation
inMedia = {
    'palmitate[s]'
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
model = setParam(model, 'lb', reactionIn(3), 0);


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

nrOfCases = 7;
fluxDistributions = zeros(length(model.rxns), nrOfCases);

for i = 1:nrOfCases
    tmpModel = model;
    switch i
        case 1
            %Fatty acid flux
            tmpModel = setParam(tmpModel, 'lb', reactionIn(1), -1);
        case 2
            %Fatty acid flux (octanoic acid)
            tmpModel = setParam(tmpModel, 'lb', reactionIn(2), -1);
        case 3
            %Canonical flux
            tmpModel = setParam(tmpModel, 'lb', reactionIn(3), -1);
        case 4
            %Complex 1 knock out
            tmpModel = setParam(tmpModel, 'lb', reactionIn(3), -1);            
            tmpModel = setParam(tmpModel, 'lb', 'HMR_6921', 0);
            tmpModel = setParam(tmpModel, 'ub', 'HMR_6921', 0);
        case 5
            %glycerol-3-phosphate shunt knock
            tmpModel = setParam(tmpModel, 'lb', reactionIn(3), -1);            
            tmpModel = setParam(tmpModel, 'lb', 'HMR_0483', 0);
            tmpModel = setParam(tmpModel, 'ub', 'HMR_0483', 0);
        case 6
            %ATP synthase KO
            tmpModel = setParam(tmpModel, 'lb', reactionIn(3), -1);            
            tmpModel = setParam(tmpModel, 'lb', 'HMR_6916', 0);
            tmpModel = setParam(tmpModel, 'ub', 'HMR_6916', 0);         
        case 7
            %oxygen transport
            tmpModel = setParam(tmpModel, 'lb', reactionIn(3), -1);            
            tmpModel = setParam(tmpModel, 'lb', 'HMR_9048', 0);
            tmpModel = setParam(tmpModel, 'ub', 'HMR_9048', 0);
        %case 8
            %ETF-ubiquinone oxidoreductase knock
            % tmpModel = setParam(tmpModel, 'lb', 'HMR_6911', 0);
            % tmpModel = setParam(tmpModel, 'ub', 'HMR_6911', 0);            
            
    end
    solution = solveLinMin(tmpModel, true);
    solution
    fluxDistributions(:, i) = solution.x;
    
end
totalSolution =  sum(abs(fluxDistributions),2);

%printAffected(model, solutionA, solutionB)

% clc
% printExchangeFluxes(model, solutionE.x)
tresh = 10^-5;
reducedModel = removeRxns(model, totalSolution<tresh, true, true, true);

model = reducedModel;
save('../reducedModel.mat', 'model');

reducedModel = rmfield(reducedModel,'rxnComps');
reducedModel = rmfield(reducedModel,'geneComps');
exportToExcelFormat(reducedModel,'../reducedModel.xls')
