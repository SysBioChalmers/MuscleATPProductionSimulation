addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;

param = loadParm('subject1', true);
param.vO2max = 1000;

settings = setSimulationSettings(2, true);
model = setupSimulation(model, param);

%Adjust stoichiometry of UCP to match ATP synthase (has a stochiometry of 3
ucp3 = findIndex(model.rxns, 'HMR_7638_m3');
model.S(:, ucp3) = 10 * model.S(:, ucp3);

%W and substrate uptake rate

reactionIn = getBounds(model, {'glycogen[s]', 'glucose[s]', 'L-lactate[s]', 'palmitate[s]'});
reactionIn(5) = getTransport(model, {'O2'}, 'sb', 's');
reactionIn(6) = getTransport(model, {'CO2'}, 'sb', 's');
reactionIn(7) = getTransport(model, {'L-lactate'}, 'sb', 'sm3');

lactateUptakeM3 = getTransport(model, {'L-lactate'}, 'sb', 'sm3');

testO2 = linspace(0,1,20);
results = zeros(length(testO2),2);

for i = 1:length(testO2)
    i
    %simulate
    model.lb(lactateUptakeM3) = -testO2(i);
    [ATPrate, maximumSolution] = runFullModel(model, settings);
    maximumSolution=maximumSolution(2,:);
    results(i,1) = molToW(1000*ATPrate(:,2));      %W max
    results(i,2) = maximumSolution(reactionIn(5)); %vO2 max
    results(i,3) = maximumSolution(ucp3);          %ucp3 flux
end


%%
clf

subplot(2,1,1)
hold all

set(gca,'xticklabel',[]) 
set(gca,'XTick',[]);

normalizedW = 100*results(:,1)/min(results(:,1));
normalizedO2 = 100*-results(:,2)/min(-results(:,2));
normalizedUCP = 100*results(:,3)/max(results(:,3));

startOfUCP = find(normalizedUCP>0);
startOfUCP = startOfUCP(1) -1;

plot(testO2, normalizedW, 'linewidth', 2)
plot(testO2, normalizedO2, 'linewidth', 2)
plot(testO2(startOfUCP) * [1 1], [100 115], 'k--', 'linewidth', 2)

ylabel('% of max')
legend({'W', 'vO2'}, 'location', 'nw')
legend boxoff
axes('Color','none','XColor','none');


subplot(2,1,2)
hold on

plot(testO2, normalizedUCP, 'linewidth', 2)
plot(testO2(startOfUCP) * [1 1], [0 100], 'k--', 'linewidth', 2)

ylabel('UCP3 (normalized flux)')
xlabel('Lactate capacity in peripheral tissue [mmol/gdw/h]')
