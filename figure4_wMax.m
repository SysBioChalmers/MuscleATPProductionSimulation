addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;
settings = [];

dwMuscle = 20*(1-0.792); %kg muscle
m1Ratio = 0.55;
vO2perDryweight = 1.38 * 1.2 * 2.8;
complex1Ratio = 33.45/76.41;
m2Efficency = 0.45;
vO2max = 11.9;
maintainance = 4;
internalWork = 2;
peripheralFA = 0;
peripheralLactateCapacity = 1.8;
FAFactor = 16.55/76.41 * 0.7;

model = setupSimulation(model, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);
settings = setSimulationSettings(20, true);

%W and substrate uptake rate
results = zeros(5,2);

reactionIn = getBounds(model, {'glycogen[s]', 'glucose[s]', 'palmitate[s]'});

%%
%W max with glycogen
maximumSolution = maximizeATP(model, settings);
wMax = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
results(1,1) = molToW(1000*wMax);
results(1,2) = maximumSolution(reactionIn(1));

%%
%W max with glucose

model.lb(reactionIn(1)) = 0;
model.lb(reactionIn(2)) = -1000;

maximumSolution = maximizeATP(model, settings);
wMax = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
results(2,1) = molToW(1000*wMax);
results(2,2) = maximumSolution(reactionIn(2));


%%
%W max with fat
peripheralFA = 1000; %we must allow periferial tissue to survive on fat
model = setupSimulation(model, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);
model.lb(reactionIn(1)) = 0;

maximumSolution = maximizeATP(model, settings);
wMax = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
results(3,1) = molToW(1000*wMax);
results(3,2) = maximumSolution(reactionIn(3));

%%
%W max with fat 2
FAFactor = 100; % we also allow unrestricted fat burning capacity
model = setupSimulation(model, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);
model.lb(reactionIn(1)) = 0;

maximumSolution = maximizeATP(model, settings);
wMax = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
results(4,1) = molToW(1000*wMax);
results(4,2) = maximumSolution(reactionIn(3));

%%
%W max with glucose uptake limitation and maximum fat

model.lb(reactionIn(1)) = 0;
model.lb(reactionIn(2)) = -100/180;

maximumSolution = maximizeATP(model, settings);
wMax = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
results(5,1) = molToW(1000*wMax);
results(5,2) = maximumSolution(reactionIn(2));


%%
%W max with "EPO"
vO2max = 1000;
model = setupSimulation(model, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor);

maximumSolution = maximizeATP(model, settings);
wMax = maximumSolution(findIndex(model.rxns, settings.primaryObjective));
results(6,1) = molToW(1000*wMax);
results(6,2) = maximumSolution(reactionIn(1));

%%
clf
hold all
glycogenStore = 0.02 * 20 * 1000/180; %2% muscle glycogen
glucoseStore = 100/180; %100 gram glycogen in liver
fatStore = 1; %some arbitary high number
totalStores = [glycogenStore; glucoseStore; fatStore; fatStore; fatStore; glycogenStore];
carbonContent = [6; 6; 16; 16; 10; 6];

wMax = results(:,1);
efficency = wMax./(carbonContent.*-results(:,2));
timeToDepletion = totalStores./-results(:,2) * 60; %minutes


timeMax = 170;

totalT = [0 timeToDepletion(1) timeToDepletion(1)+timeToDepletion(2) timeMax];
totalW = diff(totalT) .* [wMax(1) wMax(2) wMax(4)];
totalW = [0 totalW];
totalW = cumsum(totalW);


exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
%Plot epo separate
area([0 timeToDepletion(6)], [wMax(6) wMax(6)], 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none')

%Glycogen stores
area([0 totalT(2)], [wMax(1) wMax(1)], 'FaceColor', exMap(1,:), 'EdgeColor', exMap(1,:))

%Glucose stores
area([totalT(2) totalT(3)], [wMax(2) wMax(2)], 'FaceColor', exMap(2,:), 'EdgeColor', 'none')

%Glucose intake
area([totalT(3) totalT(4)], [wMax(5) wMax(5)], 'FaceColor', exMap(3,:), 'EdgeColor', 'none')

%Fat reserves
area([totalT(3) totalT(4)], [wMax(4) wMax(4)], 'FaceColor', exMap(4,:), 'EdgeColor', 'none')

%Fat reserves actual fat capacity
area([totalT(3) totalT(4)], [wMax(3) wMax(3)], 'FaceColor', exMap(5,:), 'EdgeColor', 'none')


xlabel('Minutes')
ylabel('maximum steady state W')

text(totalT(1), wMax(6), 'Epo')
text(totalT(1), wMax(1), 'Glycogen')
text(totalT(2), wMax(2), 'Glucose')
text(totalT(3), wMax(3), 'Fat capacity')
text(totalT(3), wMax(4), 'Theoretical fat capacity')
text(totalT(3), wMax(5), 'Fat + Glucose intake')
     
xlim([0 timeMax])

% color2 = [215 86 40]/256;
% plot(totalT, totalW/100, 'linewidth', 3, 'color', color2)
