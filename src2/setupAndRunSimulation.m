function [ATPrate, fullSolution] = setupAndRunSimulation(model, settings, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralFAsynth, FAFactor)
    model = addInternalConstraints(model, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, FAFactor);
    model = addTransportConstraints(model,   's', {'O2'}, [-vO2max], [0]);
    model = addTransportConstraints(model, 'sm3', {'palmitate', 'L-lactate', 'O2'}, [-peripheralFA -1000 -1],  [peripheralFAsynth 0 0]);    
    model = addTransportConstraints(model, 'sm1', {'palmitate','L-lactate'}, [-1000, -1000],  [0 0]);
    model = addTransportConstraints(model, 'sm2', {'palmitate','L-lactate'}, [-1000 0], [0 1000]);

    
    model = addSpecializedConstraints(model, 1000, maintainance, internalWork, m1Ratio);
    settings = addExchangeMedum(settings);

    settings.minRxns = {'glycogen_Exchange_', 'ventilation', 'HMR_9135_m2', 'HMR_9135_m3'};
    settings.minVal = [1, -1, -0.4, 0.01];
    [ATPrate, fullSolution] = runFullModel(model, settings);
end
