addpath('src2')
addpath('sampleData')
load('model/connectedMuscles')
model = superModel;

param = loadParm('subject1', true);
model = setupSimulation(model, param);

param.lactateBuffering = 1000; %We constrain lactate by capacity

settings = setSimulationSettings(2, false);

%
glycogenStore = 20 * 80/1000; %80 mmol/kg wet weight
subGlycogen = [param.muscleRatio 1-param.muscleRatio] * glycogenStore;
lactateCapacity = 40 * 15/1000; %40 L volume, 15 mmol concentration
glucoseStore = 0.1 * 1000/180; %100 g liver glycogen
intestinalGlucoseCapacity = 100/180; %mol/h
intestinalGlucoseCapacity = 0;

mintime = 3.7/60;
maxtime = 130/60;
timeSteps = 100;

%timePoints = 1./linspace(1/mintime, 1/maxtime, timeSteps);
timePoints1 = linspace(mintime, maxtime*0.2, timeSteps/2); %higher resolution at low times
timePoints2 = linspace(maxtime*0.21, maxtime, timeSteps/2); 
timePoints = [timePoints1 timePoints2];
timePoints = round(timePoints, 5); %to avoid numerical errors

fullSolution = zeros(length(timePoints),length(model.rxns));

internalGlycogen1 = getTransport(model, {'glycogen'}, 'sb', 'sm1');
internalGlycogen2 = getTransport(model, {'glycogen'}, 'sb', 'sm2');
internalGlucose = getTransport(model, {'glucose'}, 'sb', 's');
internalLactate = findIndex(model.rxns, 'lactateBuffering');

%W and substrate uptake rate
reactionIn = getBounds(model, {'glycogen[s]', 'glucose[s]', 'palmitate[s]'});
reactionIn(4) = getTransport(model, {'O2'}, 'sb', 's');
reactionIn(5) = getTransport(model, {'CO2'}, 'sb', 's');
reactionIn(5) = getTransport(model, {'CO2'}, 'sb', 's');
reactionIn(6) = internalLactate;


glucose=getBounds(model, {'glucose[s]'});
model.lb(glucose) = -1000;

for i = 1:length(timePoints)
    i
    %Constrain uptake rates by pool size/time
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

ATPrate = fullSolution(:,findIndex(model.rxns, settings.primaryObjective));
%saveFluxes('fluxes/duration.txt', model, fullSolution', timePoints*60, 30);
%should be run with the parsimonius flag set as true in settings

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
%Minimize O2
minO2 = zeros(length(timePoints),2);

for i = 1:length(timePoints)
    i
    %Constrain uptake rates by pool size/time
    glycogenflux1 = subGlycogen(1)/timePoints(i);
    glycogenflux2 = subGlycogen(2)/timePoints(i);
    lactateflux = lactateCapacity/timePoints(i);
    glucoseflux = glucoseStore/timePoints(i) + intestinalGlucoseCapacity;
    
    model.lb(internalGlycogen1) = -glycogenflux1;
    model.lb(internalGlycogen2) = -glycogenflux2;

    model.lb(internalGlucose) = -glucoseflux;
    model.ub(internalLactate) = lactateflux;
    minO2(i,1) = minimizeVO2(model, settings, ATPrate(i)); 
    minO2(i,2) = minimizeVO2(model, settings, ATPrate(i)*0.99);
end
%%
worldRecordData = [
    3.72	7.22
    4.75	7.02
    7.35	6.80
    12.62	6.61
    26.30	6.34
    56.43	5.91  %(allmost identical to half maraton)
    58.38	6.02
    72.42	5.88
    86.78	5.76
    122.95	5.72
    ];

estimatedWeight = 59;
steadyStateW = 324.11;

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
%normalizedW = 100*optimalW(:,1)'/max(optimalW(:,1));
normalizedW = optimalW(:,1)';


normalizedGlycogen = 6 * optimalW(:,2); %6 carbon per glycogen
normalizedGlucose = 6 * optimalW(:,3); %6 carbon per glucose
normalizedFat = 16 * optimalW(:,4); %16 carbon per fat
normalizedLactate = -3 * optimalW(:,7); %3 carbon per lactate

totalCarbon = normalizedGlycogen + normalizedGlucose + normalizedFat;

Yvalues = [normalizedW .* normalizedGlycogen'./totalCarbon';
           normalizedW .* normalizedGlucose'./totalCarbon';
           normalizedW .* normalizedFat'./totalCarbon'];

Yvalues = [Yvalues(1,:); Yvalues];       
       
Yvalues(1,:) =  (1-(normalizedLactate./normalizedGlycogen))' .* Yvalues(1,:);
Yvalues(2,:) =  (normalizedLactate./normalizedGlycogen)' .* Yvalues(2,:);

tmp = Yvalues(2,:);
Yvalues(2,:) = Yvalues(1,:);
Yvalues(1,:) = tmp;


h = area(timePoints*60, Yvalues', 'FaceColor', exMap(1,:), 'EdgeColor', 'none');
h(1).FaceColor  = color2;
h(3).FaceColor  = exMap(2,:);
h(4).FaceColor  = exMap(3,:);



plot(timePoints*60,normalizedW,'linewidth',2)

plot([0 max(timePoints)*60], steadyStateW * [1 1], 'k-');
%plot(timePoints*60, 100*steadyStateLactate/max(steadyStateLactate), 'k--', 'linewidth', 2);
%plot(worldRecordData(:,1), 100*worldRecordData(:,2), 'kx', 'linewidth', 2);

xlim(60*[min(timePoints) max(timePoints)])
ylim([0 inf])
xlabel('minutes')
ylabel('work rate [W]')


legend('glycogen (to lactate)', 'glycogen', 'glucose', 'fat',  'watt', 'location', 'ne')
legend boxoff
xlabel('minutes')

figure()
hold all
vO2Optimal = -optimalW(:,5)/max(-optimalW(:,5));
vO2NearOptimal = -minO2./max(-optimalW(:,5));
divergence = optimalW(:,3)<-0.001;
divergence = find(divergence)-1;
divergenceT = timePoints(divergence(1)) * 60;
area(timePoints*60, 100* vO2NearOptimal(:,2)', 'FaceColor', exMap(2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off')

plot(timePoints*60, 100 * vO2Optimal,'-', 'color', exMap(1,:), 'linewidth',2)
plot(timePoints*60, 99 * vO2Optimal,'--', 'color', exMap(1,:),  'linewidth',1)
plot(timePoints*60, 100 * vO2NearOptimal(:,2),'-', 'color', exMap(2,:),  'linewidth',2)    
plot(divergenceT * [1 1], [0 100], 'k--')

timeInterval = and(timePoints*60<115, timePoints*60>14);
mean(1-vO2NearOptimal(timeInterval,2)./vO2Optimal(timeInterval))
min(1-vO2NearOptimal(timeInterval,2)./vO2Optimal(timeInterval))
max(1-vO2NearOptimal(timeInterval,2)./vO2Optimal(timeInterval))

xlim(60*[min(timePoints) max(timePoints)])
ylim([70 100])
xlabel('minutes')
ylabel('%of max')

legend('100% Wmax', '100% Wmax * (.99)', '99% Wmax', 'location', 'ne')
legend boxoff

figure()

subplot(5,1,2:5)
hold all

wattEstimates = interp1(timePoints,optimalW(:,1),worldRecordData(:,1)/60);
wattPerKg = wattEstimates/estimatedWeight;

for i = 1:length(worldRecordData(:,2))
%   text(wattPerKg(i), worldRecordData(i,2)+0.05, sprintf('%2.1f min (%2.1fKm)', worldRecordData(i,1), distance(i)));
    text(wattPerKg(i), worldRecordData(i,2)+0.05, sprintf('%2.1f', worldRecordData(i,1)));   
end
[r, p] = corr(wattPerKg, worldRecordData(:,2))

xvals = linspace(4, 8);
mdl = fitlm(wattPerKg, worldRecordData(:,2),'Intercept',false);
plot(xvals, mdl.Coefficients.Estimate(1) * xvals, '-', 'color', [0.5 0.5 0.5]); 
mdl

mdl = fitlm(wattPerKg, worldRecordData(:,2));
plot(xvals, mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * xvals, '-', 'color', exMap(1,:), 'linewidth', 2);
mdl

scatter(wattPerKg, worldRecordData(:,2), 'filled', 'markerfacecolor', color2)

xlim([4 8])
ylim([5.5 7.5])
xlabel('predicted specific work rate [W/kg]')
ylabel('world record average speed [m/s]')

subplot(5,1,1)
hold all
area([wattPerKg; 0], [distance; 1000], 'LineStyle', 'none', 'FaceColor', [0.6 0.6 0.6])
scatter(wattPerKg, distance, 'filled')
xlabel('duration [min]')
xlim([4 8])
minMax = round([min(distance) max(distance)],1);
ylim(minMax)
yticks(minMax)
%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})

ax1 = gca; % current axes
ax1.XColor = 'none';
