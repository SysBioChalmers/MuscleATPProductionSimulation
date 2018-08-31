function model = setupSimulation(model, maintainance, internalWork, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, vO2max, peripheralFA, peripheralLactateCapacity, FAFactor, type1tresh, type2target)
    model = addInternalConstraints(model, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, FAFactor, type1tresh, type2target);
    model = addTransportConstraints(model,   's', {'O2'}, [-vO2max], [0]);
    model = addTransportConstraints(model, 'sm3', {'palmitate', 'L-lactate'}, [-peripheralFA -peripheralLactateCapacity],  [0 0]);    
    model = addTransportConstraints(model, 'sm1', {'palmitate','L-lactate'}, [-1000, -1000],  [0 0]);   
    model = addTransportConstraints(model, 'sm2', {'palmitate','L-lactate'}, [-1000 0], [0 1000]);
    model = addMaintenanceConstraints(model, 1000, maintainance, internalWork);
    
    %Allow production of stearate
    mTransp = getTransport(model, {'stearate'}, 's', 'sb');
    model.ub(mTransp) = 1000;
    mTransp = getTransport(model, {'stearate'}, 'sm3', 'sb');
    model.ub(mTransp) = 1000;
    
    model = addExchangeMedum(model);

end

