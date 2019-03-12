addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;

param = loadParm('subject1', true);
param.lactateBuffering = 0;

settings = setSimulationSettings(2, false);
model = setupSimulation(model, param);

%W and substrate uptake rate
results = zeros(5,2);

reactionIn = getBounds(model, {'glycogen[s]', 'glucose[s]', 'L-lactate[s]', 'palmitate[s]'});
reactionIn(5) = getTransport(model, {'O2'}, 'sb', 's');
reactionIn(6) = getTransport(model, {'CO2'}, 'sb', 's');
reactionIn(7) = getTransport(model, {'L-lactate'}, 'sb', 'sm2');

lactateUptakeM2 = getTransport(model, {'L-lactate'}, 'sb', 'sm2');

model.lb(reactionIn(1:4)) = 0;

modelRef = model;
paramRef = param;
caseLabels = {
    'glycogen'
    'glucose'
    'lactate'
    'fat (unconstrained)'
    'fat'
    'glycogen (unconstrained O2)'};

for i = 1:length(caseLabels)
    switch i
    case 1
        %W max with glycogen
        model.lb(reactionIn(1)) = -1000;
    case 2
        %W max with glucose
        model.lb(reactionIn(2)) = -1000;   
    case 3
        %W max with lactate
        param.peripheralLactateCapacity = 1000; %we must allow periferial tissue to survive on lactate        
        model = setupSimulation(model, param);
        model.lb(lactateUptakeM2) = -1000; %we must allow m2 to uptake laktate
        model.lb(reactionIn(1:4)) = 0;
        model.lb(reactionIn(3)) = -1000;   
    case 4
        %W max with fat
        param.peripheralFA = 1000; %we must allow periferial tissue to survive on fat
        param.HMR_6911 = 1000; % we allow unrestricted fat burning capacity
        model = setupSimulation(model, param);
        model.lb(reactionIn(1:4)) = 0;
        model.lb(reactionIn(4)) = -1000;   
    case 5
        %W max with fat
        param.peripheralFA = 1000; %we must allow periferial tissue to survive on fat
        model = setupSimulation(model, param);
        model.lb(reactionIn(1:4)) = 0;
        model.lb(reactionIn(4)) = -1000;   
    case 6
        param.vO2max = 1000;
        model = setupSimulation(model, param);
        model.lb(reactionIn(1)) = -1000;
    end
    
    %simulate
    [ATPrate, maximumSolution] = runFullModel(model, settings);
    maximumSolution=maximumSolution(2,:);
    results(i,1) = molToW(1000*ATPrate(:,2));      %W max
    results(i,2) = maximumSolution(reactionIn(5)); %vO2 max
    
    %reset parameters
    model = modelRef;
    param = paramRef;    
end


%%
clf
hold all
wMax = results(:,1);

exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
     
for i = 1:(length(caseLabels)-1)
    barN = length(caseLabels)-i+1;
    barh(barN, wMax(i), 'FaceColor', exMap(i,:), 'EdgeColor', exMap(i,:))
    text(1,barN, sprintf('%2.0f%%',100*wMax(i)/wMax(1)));
end

%Plot epo separate
i = length(caseLabels);
barh(1, wMax(i), 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none')
text(1,1, sprintf('%2.0f%%',100*wMax(i)/wMax(1)));
text(100,1, sprintf('(vO2 %2.0f%%)',100*results(i,2)/results(1,2)));

yticks(1:length(caseLabels))
yticklabels(flipud(caseLabels))

ylim([0.5 length(caseLabels) + 0.5])
xlim([0 450])
xlabel('maximum steady state W')
plot(wMax(1)*[1 1], [0.5 length(caseLabels) + 0.5], 'k--')

dO2 = molToMl(-(results(5,2)-results(1,2)));
dW = results(5,1)-results(1,1);
dO2dw = dO2/dW

xlim([0 400])

%Fat capacity:
-results(4,2)*256.43/60
-results(3,2)*256.43/60




