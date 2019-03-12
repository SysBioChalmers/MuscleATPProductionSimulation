addpath('../src2')
addpath('../sampleData')
load('../model/connectedMuscles')
model = superModel;

param = loadParm('subject1', true);

settings = setSimulationSettings(2, false);
model = setupSimulation(model, param);

%W and substrate uptake rate
reactionIn = getBounds(model, {'L-lactate[s]'});
lactateInM2 = getTransport(model, {'L-lactate'}, 'sb', 'sm2');
model.lb(lactateInM2) = -1000;

lactateFlux = linspace(0,6, 10);
results = zeros(length(lactateFlux),1);

for i = 1:length(lactateFlux)
    model.lb(reactionIn) = -lactateFlux(i);
    model.ub(reactionIn) = -lactateFlux(i);

    [ATPrate, maximumSolution] = runFullModel(model, settings);
    if length(maximumSolution) == 1
        ATPrate = 0;
    else
        ATPrate = molToW(1000*ATPrate(2));
    end
    
    results(i,1) = ATPrate;      %W max
end

%%
hold all
aproxlactateSS = 2;
aproxWmax = interp1q(lactateFlux, results, aproxlactateSS);
SSwMax = aproxWmax/max(results);

plot(lactateFlux, results, 'linewidth', 2)
% plot(aproxlactateSS * [0 1 1], aproxWmax * [1 1 0], 'k--')
% text(aproxlactateSS, aproxWmax, sprintf('%2.2f%% of W max ', 100*SSwMax));

xlabel('Lactate flux [mol/gdw/h]')
ylabel('W max')
xlim([0 4])
ylim([0 inf])
