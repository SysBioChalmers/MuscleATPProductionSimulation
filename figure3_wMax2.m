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
reactionIn(6) = getTransport(model, {'L-lactate'}, 'sb', 'sm3');
reactionIn(7) = getTransport(model, {'L-lactate'}, 'sb', 's');

%
glycogenStore = 20 * 80/1000; %80 mmol/kg wet weight
subGlycogen = [muscleRatio 1-muscleRatio] * glycogenStore;
lactateCapacity = 37 * 15/1000;
glucoseStore = 0.1 * 1000/180; %100 g liver glycogen
intestinalGlucoseCapacity = 100/180; %mol/h

mintime = 3.716/60;
maxtime = 130/60;
timeSteps = 100;

%timePoints = 1./linspace(1/mintime, 1/maxtime, timeSteps);
timePoints1 = linspace(mintime, maxtime*0.2, timeSteps/2); 
timePoints2 = linspace(maxtime*0.25, maxtime, timeSteps/2); 
timePoints = [timePoints1 timePoints2];


fullSolution = zeros(length(timePoints),length(model.rxns));

internalGlycogen1 = getTransport(model, {'glycogen'}, 'sb', 'sm1');
internalGlycogen2 = getTransport(model, {'glycogen'}, 'sb', 'sm2');
internalGlucose = getTransport(model, {'glucose'}, 'sb', 's');
internalLactate = getTransport(model, {'L-lactate'}, 'sb', 's');

lactateFlux=getBounds(model, {'L-lactate[s]'});
model.ub(lactateFlux) = 1000;

palmitate=getBounds(model, {'palmitate[s]'});
model.lb(palmitate) = -1000;

glucose=getBounds(model, {'glucose[s]'});
model.lb(glucose) = -1000;

%oxygen delivery to muscle may be limiting
model = limitOxygenDelivery(model, internalO2);

for i = 1:length(timePoints)
    i
    glycogenflux1 = subGlycogen(1)/timePoints(i);
    glycogenflux2 = subGlycogen(2)/timePoints(i);
    lactateflux = lactateCapacity/timePoints(i);
    glucoseflux = glucoseStore/timePoints(i) + intestinalGlucoseCapacity;
    
    model.lb(internalGlycogen1) = -glycogenflux1;
    model.lb(internalGlycogen2) = -glycogenflux2;

    model.lb(internalGlucose) = -glucoseflux;
    model.ub(internalLactate) = lactateflux;
       
    
    %maximumSolution = maximizeATP(model, settings);
    [ATPrate, maximumSolution] = runFullModel(model, settings);
    fullSolution(i,:) = maximumSolution(2,:)';

end

%Remove the oxygen tradeoff psuedo metabolite after simulation since it 
%interferes with data retrival
model.S(findIndex(model.mets, 'oxygenTradeof'),:) = 0;

ATPrate = fullSolution(:,findIndex(model.rxns, settings.primaryObjective));

optimalW = zeros(length(timePoints),length(reactionIn)+1);
optimalW(:,1) = molToW(1000*ATPrate);
optimalW(:,2:end) = fullSolution(:, reactionIn);

%     optimalW(i,2) = maximumSolution(internalGlycogen1);
%     optimalW(i,3) = maximumSolution(internalGlycogen2);
%     optimalW(i,4) = maximumSolution(internalGlycogen3);
%     optimalW(i,5) = maximumSolution(internalLactate);

steadyStateLactate = convertVlacToConc(-optimalW(:,7), 1.4, 9);

fprintf('Time [h]\tW\tglycogen\tglucose	palmitate\tO2\tCO2\tLactate flux\tLactate concentration\n')

for i = 1:length(timePoints)
    fprintf('%2.2f', timePoints(i)) 
    for j = 1:size(optimalW,2)
        fprintf('\t%2.2f', optimalW(i,j)); 
    end
    fprintf('\t%2.2f\n', steadyStateLactate(i)); 
end
%%
worldRecordData = [
    3.72	7.22
    4.75	7.02
    7.35	6.80
    12.62	6.61
    26.30	6.34
    56.43	5.91  %allmost identical to half maraton
    58.38	6.02
    72.42	5.88
    86.78	5.76
    122.95	5.72
    ];

estimatedWeight = 59;

distance = worldRecordData(:,1).*worldRecordData(:,2)*60/1000;

%worldRecordData(:,2) = worldRecordData(:,2)/max(worldRecordData(:,2));

close all
exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
color2 = [215 86 40]/256;


figure()
hold all
normalizedW = 100*optimalW(:,1)'/max(optimalW(:,1));

normalizedGlycogen = 6 * optimalW(:,2); %6 carbon per glycogen
normalizedGlucose = 6 * optimalW(:,3); %6 carbon per glucose
normalizedFat = 16 * optimalW(:,4); %16 carbon per fat
normalizedLactate = -3 * optimalW(:,8); %3 carbon per lactate

totalCarbon = normalizedGlycogen + normalizedGlucose + normalizedFat;

Yvalues = [normalizedW .* normalizedGlycogen'./totalCarbon';
           normalizedW .* normalizedGlucose'./totalCarbon';
           normalizedW .* normalizedFat'./totalCarbon'];

Yvalues = [Yvalues(1,:); Yvalues];       
       
Yvalues(1,:) =  (1-(normalizedLactate./normalizedGlycogen))' .* Yvalues(1,:);
Yvalues(2,:) =  (normalizedLactate./normalizedGlycogen)' .* Yvalues(2,:);



h = area(timePoints*60, Yvalues', 'FaceColor', exMap(1,:), 'EdgeColor', 'none');
h(2).FaceColor  = color2;
h(3).FaceColor  = exMap(2,:);
h(4).FaceColor  = exMap(3,:);



plot(timePoints*60,100*optimalW(:,1)'/max(optimalW(:,1)),'linewidth',2)
plot(timePoints*60,100* -optimalW(:,5)/max(-optimalW(:,5)),'k-', 'linewidth',2)

%plot(timePoints*60, 100*steadyStateLactate/max(steadyStateLactate), 'k--', 'linewidth', 2);
%plot(worldRecordData(:,1), 100*worldRecordData(:,2), 'kx', 'linewidth', 2);

xlim(60*[min(timePoints) max(timePoints)])
ylim([0 inf])
xlabel('minutes')
ylabel('%of max')


legend('glycogen', 'lactate (accumulated)', 'glucose', 'fat',  'watt', 'vO2', 'World records', 'location', 'se')
xlabel('minutes')
figure()
hold all


wattEstimates = interp1(timePoints,optimalW(:,1),worldRecordData(:,1)/60);
wattPerKg = wattEstimates/estimatedWeight;

% mdl1 = fitlm(worldRecordData(:,2), speedEstimates);
% km = mdl1.Coefficients.Estimate;
% 
% plotInterval = [20 27];

% corr( worldRecordData(:,2), speedEstimates)
% plot(plotInterval, km(1) + km(2) * plotInterval, 'k--', 'linewidth', 2);
% plot(plotInterval, plotInterval, 'linewidth', 2, 'color', 0.5*[1 1 1]);
% scatter(worldRecordData(:,2), speedEstimates, 'filled')
% 
% xlabel('World record average speed [km/h]')
% ylabel('Predicted optimal speed [km/h]')

%scatter(worldRecordData(:,1), worldRecordData(:,2), 'filled')

for i = 1:length(worldRecordData(:,2))
   text(worldRecordData(i,2)+0.05, wattPerKg(i), sprintf('%2.1f min (%2.1fKm)', worldRecordData(i,1), distance(i)));
end

scatter(worldRecordData(:,2), wattPerKg, 'filled')
%set(gca, 'XScale', 'log')
xlim([5.5 7.5])
ylim([4 8])
xlabel('World record average speed [m/s]')
ylabel('Predicted optimal specific power [W/kg]')
legend boxoff  

%selectedTimes = [1 3 5 7 10];

%scatter(wattEstimates(selectedTimes), worldRecordData(selectedTimes,2), 'filled')
%xlim(plotInterval)
%ylim([14 24])

% figure()
% compareWithSampleData(model, 'subject1', ATPrate, fullSolution)
% figure()
% plotMetaboliteList = {'O2', 'glycogen', 'L-lactate', 'CO2', 'glucose', 'palmitate', 'stearate'};
% plotFullSolutionInternalBlood(model, ATPrate, fullSolution, plotMetaboliteList, 'sb', {'sm1', 'sm2', 'sm3'});
% %figure()
% %simData = plotFullSolutionAandB(model, ATPrate, fullSolution, 'sb', {'sm1', 'sm2'});
% figure();
% plotFiberDistribution(model, ATPrate', fullSolution, muscleRatio)
