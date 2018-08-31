addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;

dwMuscle = 20*(1-0.792); %kg muscle
muscleRatio = 0.55; %type1
vO2perDryweight = 1.38 * 2.6; %Specific activity Complex IV
complex1Ratio = 33.45/76.41;
m2Efficency = 0.5;
vO2max = 1000;
internalO2 = 10;
maintainance = 5.5;
internalWork = 2;
peripheralFA = 0.004;
peripheralLactateCapacity = 2;
FAFactor = 0.8*16.55/76.41;
type1tresh = 9; %type 1 is activated first, 10 mol ATP corresponds to 20% of Wmax
type2target = 0.55; %type 2 dominates at the highest W rates

model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);
settings = setSimulationSettings(2, false);


%W and substrate uptake rate
results = zeros(5,3);

reactionIn = getBounds(model, {'glycogen[s]', 'glucose[s]', 'palmitate[s]'});
reactionIn(4) = getTransport(model, {'O2'}, 'sb', 's');
reactionIn(5) = getTransport(model, {'CO2'}, 'sb', 's');
reactionIn(6) = getTransport(model, {'L-lactate'}, 'sb', 'sm2');

%oxygen delivery to muscle may be limiting
model = limitOxygenDelivery(model, internalO2);


%%
%W max with glycogen
[ATPrate, maximumSolution] = runFullModel(model, settings);
maximumSolution=maximumSolution(2,:);
results(1,1) = molToW(1000*ATPrate(:,2));
%results(1,2) = maximumSolution(reactionIn(1));
results(1,3) = maximumSolution(reactionIn(4));

%%
%W max with glucose

model.lb(reactionIn(1)) = 0;
model.lb(reactionIn(2)) = -1000;

[ATPrate, maximumSolution] = runFullModel(model, settings);
maximumSolution=maximumSolution(2,:);
results(2,1) = molToW(1000*ATPrate(:,2));
% results(2,2) = maximumSolution(reactionIn(2));
% results(2,3) = maximumSolution(reactionIn(4));


%%
%W max with fat
model.S(end,:) = 0;
peripheralFA = 1000; %we must allow periferial tissue to survive on fat
model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);
model.lb(reactionIn(1)) = 0;
model.lb(reactionIn(2)) = 0;
model = limitOxygenDelivery(model, internalO2);

[ATPrate, maximumSolution] = runFullModel(model, settings);
maximumSolution=maximumSolution(2,:);
results(3,1) = molToW(1000*ATPrate(:,2));
results(3,2) = maximumSolution(reactionIn(3));
% results(3,3) = maximumSolution(reactionIn(4));

%%
%W max with fat 2
peripheralFA = 1000; %we must allow periferial tissue to survive on fat
FAFactor = 10; % we also allow unrestricted fat burning capacity
model.S(end,:) = 0; %remove oxygen constraint
model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);
model.lb(reactionIn(1)) = 0;
model.lb(reactionIn(2)) = 0;
model = limitOxygenDelivery(model, internalO2);

[ATPrate, maximumSolution] = runFullModel(model, settings);
maximumSolution=maximumSolution(2,:);
results(4,1) = molToW(1000*ATPrate(:,2));
results(4,2) = maximumSolution(reactionIn(3));
% results(4,3) = maximumSolution(reactionIn(4));

%%
%W max with "EPO"
vO2max = 1000;
model.S(end,:) = 0;
model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);


[ATPrate, maximumSolution] = runFullModel(model, settings);
maximumSolution=maximumSolution(2,:);
results(5,1) = molToW(1000*ATPrate(:,2));
results(5,2) = maximumSolution(reactionIn(1));
results(5,3) = maximumSolution(reactionIn(4));

%%
clf
hold all
wMax = results(:,1);

exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
%Glycogen
barh(5, wMax(1), 'FaceColor', exMap(1,:), 'EdgeColor', exMap(1,:))
text(1,5, sprintf('%2.0f%%',100*wMax(1)/wMax(1)));

%Glucose
barh(4, wMax(2), 'FaceColor', exMap(2,:), 'EdgeColor', 'none')
text(1,4, sprintf('%2.0f%%',100*wMax(2)/wMax(1)));

%Fat
barh(3, wMax(4), 'FaceColor', exMap(3,:), 'EdgeColor', 'none')
text(1,3, sprintf('%2.0f%%',100*wMax(4)/wMax(1)));

%Fat, actual fat capacity
barh(2, wMax(3), 'FaceColor', exMap(4,:), 'EdgeColor', 'none')
text(1,2, sprintf('%2.0f%%',100*wMax(3)/wMax(1)));

%Plot epo separate
barh(1, wMax(5), 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none')
text(1,1, sprintf('%2.0f%%',100*wMax(5)/wMax(1)));
text(100,1, sprintf('(vO2 %2.0f%%)',100*results(5,3)/results(1,3)));


yticks(1:5)
yticklabels({'Unconstrained O2', 'Fat capacity', 'Theoretical fat capacity', 'Glucose', 'Glycogen'})
ylim([0.5 5.5])
xlim([0 450])
xlabel('maximum steady state W')
plot(wMax(1)*[1 1], [0.5 5.5], 'k--')

dO2 = molToMl(-(results(5,3)-results(1,3)));
dW = results(5,1)-results(1,1);
dO2dw = dO2/dW

xlim([0 400])

%Fat capacity:
-results(4,2)*256.43/60
-results(3,2)*256.43/60




