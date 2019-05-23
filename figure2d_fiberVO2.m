close all
hold all
load('model/reducedModel')
addpath('src1')
addpath('model/makeMultiTissueModel')
saturation = 0.5;

color3 = [237 197 188]/256;
color2 = [215 86 40]/256;
color1 = [93 155 211]/256;

model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = mapProteomToRxns(model, 'data/RxnAndProtein.txt');

%We do not trust proteomics for 6-PHOSPHOFRUCTOKINASE
model.proteinMass(findIndex(model.rxns,  'HMR_4379')) = -1;

[model, constraindRxns] = addProteinConstrant(model, saturation);
objectiveFunction = {'human_ATPMaintainance'};
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);
minimalMedia = {
    'O2[s]'
    'octanoic acid[s]'
    'glycogen[s]'
    'H2O[s]'
    };

minimalFlux= [-1000
              0
              0
              -1000];

reactionNumbers = getBounds(model, minimalMedia);
model = setParam(model, 'lb', reactionNumbers, minimalFlux);


rxnsToAdd = createRXNStuct(model, 'glutIn', 'aspartate[c] => glutamate[c]', 0, 0, 'Glutamate as substrate');
model=addRxns(model,rxnsToAdd,3,'sb',false);

rxnsToAdd = createRXNStuct(model, 'sucIn', 'fumarate[m] => succinate[m]', 0, 0, 'Succinate as substrate');
model=addRxns(model,rxnsToAdd,3,'sb',false);

%simulate VO2
simulationData = zeros(5,2);

for i = 1:length(simulationData)
    tmpModel = model;
    
    switch i
    case 1
        %VO2 fat
        tmpModel = setParam(tmpModel, 'lb', reactionNumbers(2), -1000);
    case 2
        %VO2 glutamate
        tmpModel = setParam(tmpModel, 'ub', 'glutIn', 1000);
    case 3
        %VO2 glutamate + succinate
        tmpModel = setParam(tmpModel, 'ub', 'glutIn', 1000);
        tmpModel = setParam(tmpModel, 'ub', 'sucIn', 1000);
    case 4
        %VO2 glutamate + succinate + fat
        tmpModel = setParam(tmpModel, 'lb', reactionNumbers(2), -1000);
        tmpModel = setParam(tmpModel, 'ub', 'glutIn', 1000);
        tmpModel = setParam(tmpModel, 'ub', 'sucIn', 1000);
    case 5
        %VO2 glutamate + succinate (uncoupling)
        tmpModel = setParam(tmpModel, 'ub', 'glutIn', 1000);
        tmpModel = setParam(tmpModel, 'ub', 'sucIn', 1000);
        tmpModel.ub(findIndex(tmpModel.rxns, 'HMR_6921')) = 1000;%CI
        %tmpModel.lb(findIndex(tmpModel.rxns, 'HMR_4652')) = -1000;%CII
        tmpModel.ub(findIndex(tmpModel.rxns, 'HMR_6918')) = 1000;%CIII
        tmpModel.ub(findIndex(tmpModel.rxns, 'HMR_6914')) = 1000;%CIV
        tmpModel.ub(findIndex(tmpModel.rxns, 'HMR_6916')) = 1000;%CV
    end
        
    solution = solveLinMin(tmpModel);
    simulationData(i, 1) = -solution.x(reactionNumbers(1));

    tmpModel.lb(findIndex(tmpModel.rxns, 'HMR_0483')) = 0; %block GLY-PHOS    
    solution = solveLinMin(tmpModel);
    simulationData(i, 2) = -solution.x(reactionNumbers(1));
    
end


%%   
factor = 1/(1-0.792) * 60*60 * 10^-12 * 10^3 * 10^3; %per hour per gdw 

experimentalData = factor*[
    16.55   1.41 % Fat
    33.45	2.11 % Complex I
%    26.76	1.37 % Complex I + pyruvate
    76.41	4.58 % Complex I+II
    69.88   3.92 % Complex I+II + fat
    105.99	5.63 % FCCP
    ];

     
exMap = [color1;
         color2;
         color3
         ];

     

     
hold all     
points = 1:2:(2*length(simulationData) -1);

plotData = simulationData(:,[2 1]);
plotData(:,2) = plotData(:,2) - plotData(:,1);


h = bar(points, plotData, 0.35, 'stacked', 'FaceColor', exMap(2,:), 'EdgeColor', 'none');
h(2).FaceColor = exMap(3,:);

bar(points+1, experimentalData(:,1), 0.35, 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
errorbar(points+1, experimentalData(:,1), experimentalData(:,2),'k.')
%bar([2 4], [Complex1O2 Complex2O2], 0.5, 'FaceColor', exMap(2,:), 'EdgeColor', 'none')

set(gca, 'XTick', points+0.5)
set(gca, 'XTickLabel', {'F', 'CI', 'CI+II', 'CI+II+F', 'CI+II (u)'})
ylabel('vO2 [mmol/gdw/h]')
ylim([0, 2.5])
xlim([0.5 10.5])
xtickangle(45)
set(findall(gcf,'-property','FontSize'),'FontSize',13)
legend({'model', 'model (complex I bypass)', 'data'}, 'location', 'NW')
legend boxoff