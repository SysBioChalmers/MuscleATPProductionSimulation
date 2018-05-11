function [ATPrate, fullSolution] = setupAndRunSimulation(model, settings, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor)
    model = addInternalConstraints(model, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, FAFactor);
    model = addTransportConstraints(model,   's', {'O2'}, [-vO2max], [0]);
    model = addTransportConstraints(model, 'sm3', {'palmitate', 'L-lactate'}, [-peripheralFA -peripheralLactateCapacity],  [0 0]);    
    model = addTransportConstraints(model, 'sm1', {'palmitate','L-lactate'}, [-1000, -1000],  [0 0]);
    model = addTransportConstraints(model, 'sm2', {'palmitate','L-lactate'}, [-1000 0], [0 1000]);
     
    model = addSpecializedConstraints(model, 1000, maintainance, internalWork, m1Ratio);
    settings = addExchangeMedum(settings);

    %Allow production of decanoic acid consumptuon
    mTransp = getTransport(model, {'stearate'}, 's', 'sb');
    model.ub(mTransp) = 1000;
    mTransp = getTransport(model, {'stearate'}, 'sm3', 'sb');
    model.ub(mTransp) = 1000;


    settings.minRxns = {'glycogen_Exchange_', 'O2_Exchange_b', 'HMR_9135_m2', 'HMR_9135_m3'};
    settings.minVal = [1, 1, -0.4, 0.3];
    [ATPrate, fullSolution] = runFullModel(model, settings);

end

