close all
hold all
load('model/reducedModel')
addpath('src1')
saturation = 0.5;

model = mapDataToRxns(model, 'data/RxnAndSA.txt');
model = mapProteomToRxns(model, 'data/RxnAndProtein.txt');

%We do not trust proteomics for 6-PHOSPHOFRUCTOKINASE
model.proteinMass(findIndex(model.rxns,  'HMR_4379')) = -1;

%We do not trust proteomics for TPI
model.proteinMass(findIndex(model.rxns,  'HMR_4391')) = -1;



[model, constraindRxns] = addProteinConstrant(model, saturation);

minimalMedia = {
    'O2[s]'
    'octanoic acid[s]'
    'glycogen[s]'
    'H2O[s]'
    };

plotExchange = {
    'O2[s]'
    'octanoic acid[s]'    
    'glycogen[s]'
    'L-lactate[s]'
    };

minimalFlux= [-1000
              -1000
              -1000
              -1000];


reactionNumbers = getBounds(model, minimalMedia);

model = setParam(model, 'lb', reactionNumbers, minimalFlux);

objectiveFunction = {'human_ATPMaintainance'};
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);
model = setParam(model, 'obj', objectiveFunction, 1);
solution = solveLinMin(model);
maxATP = -solution.f;
atpRates = linspace(0, maxATP, 1000);

model = setParam(model, 'obj', reactionNumbers(3), 1);
fullSolution = runChemostatExperiment(model, atpRates, objectiveFunction);
%saveFluxes('fluxes/fiber.txt', model, fullSolution', atpRates, 30);


oxygenFlux = fullSolution(:, findIndex(model.rxns, 'HMR_9048'));
DeltaO2vsDeltaATP = (diff(oxygenFlux)./diff(atpRates'));
DeltaO2vsDeltaATP(20)/DeltaO2vsDeltaATP(2)
breakPoints = find(abs(diff(DeltaO2vsDeltaATP))>0.01);
breakPoints = breakPoints([1 2 6]) + 1;
breakPoints = [1; breakPoints; length(oxygenFlux)];



%%
hold all
fillColors = [232 208 190
              190 209 223    
              199 227 187
              233 233 233
            ]/256;

for i = 1:(length(breakPoints)-1)
    curInterval = breakPoints(i):breakPoints(i+1);
    area(atpRates(curInterval),-oxygenFlux(curInterval), 'FaceColor', fillColors(i,:), 'EdgeColor', 'none', 'HandleVisibility','off')
end


cmolSoluiton = fullSolution;
cmolSoluiton(:,reactionNumbers(2)) = cmolSoluiton(:,reactionNumbers(2))*2;
plotFullSolution(model, atpRates, cmolSoluiton, plotExchange);

% for i = 1:length(breakPoints)
%     plot(atpRates(breakPoints(i)) * [0 1], -oxygenFlux(breakPoints(i))*[1 1], 'k--', 'HandleVisibility','off')
%     plot(atpRates(breakPoints(i)) * [1 1], -oxygenFlux(breakPoints(i))*[0 1], 'k--', 'HandleVisibility','off')
% end
I = breakPoints(2):breakPoints(3);
slope = polyfit(atpRates(I), -oxygenFlux(I)',1);
I = breakPoints(2):breakPoints(4);
plot(atpRates(I), atpRates(I)*slope(1) + slope(2), 'k--', 'HandleVisibility','off')



xlim([0 8])
ylim([0 2])


