load('model/connectedMuscles')
addpath('src2')
addpath('sampleData')
model = superModel;


%Random sampling settings
pertubationFactor = 1.001; %.1 difference


%Parameterset
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

%Make reference condition
model = setupSimulation(model, maintainance, internalWork, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target);
settings = setSimulationSettings(2, false);

%Reactions to plot
readouts = getTransport(model, {'O2', 'CO2', 'glycogen'}, 'sb', 's');

%apply oxygen constraint
model = limitOxygenDelivery(model, internalO2);

%Make reference condition
[referenceATP, referenceSolution] = runFullModel(model, settings);

%Store parameters
parameterList = zeros(14,1);
parameterList(1) = maintainance;
parameterList(2) = internalWork;
parameterList(3) = dwMuscle;
parameterList(4) = muscleRatio;
parameterList(5) = vO2perDryweight;
parameterList(6) = m2Efficency;
parameterList(7) = complex1Ratio;
parameterList(8) = vO2max;
parameterList(9) = peripheralFA;
parameterList(10) = peripheralLactateCapacity;
parameterList(11) = FAFactor;
parameterList(12) = type1tresh;
parameterList(13) = type2target;
parameterList(14) = internalO2;

parameterNames = {'maintainance'
                  'internalWork'
                  'dwMuscle'
                  'm1Ratio'
                  'vO2perDryweight'
                  'm2Efficency'
                  'complex1Ratio'
                  'vO2max'
                  'peripheralFA'
                  'peripheralLactateCapacity'
                  'FAFactor'
                  'Type 1 tresh'
                  'Type 2 target'
                  'internal O2'};


%Result vectors
results = zeros(length(parameterList),4);


%Reference result
resultReference = zeros(1,4);
resultReference(1) = referenceATP(end);
resultReference(2:4) = referenceSolution(end, readouts);
resultReference(5) = resultReference(3)/-resultReference(2);

for i = 1:length(parameterList)
    perturb = parameterList;
    perturb(i) = perturb(i)*pertubationFactor;
    model = setupSimulation(superModel, perturb(1), perturb(2), perturb(3),  perturb(4), perturb(5), perturb(6), perturb(7), perturb(8), perturb(9), perturb(10), perturb(11), perturb(12), perturb(13));
    model = limitOxygenDelivery(model, perturb(14));
    [ATPrate, fullSolution] = runFullModel(model, settings);
    results(i,1) = ATPrate(end);
    results(i,2:4) = fullSolution(end, readouts);  
    results(i,5) = results(i,3)/-results(i,2);
end

%%
%Plot sampling result
fprintf('parameterName\tWatt\tvO2\tvCo2\tvglycogen\tRQ\n')

for i = 1:length(parameterList)
    percentCange = resultReference./results(i,:);
    E = -(1-percentCange)./(1-pertubationFactor);
    fprintf('%s\t%2.5f\t%2.5f\t%2.5f\t%2.5f\t%2.5f\n', parameterNames{i}, E(1), E(2), E(3), E(4), E(5)); 
end

