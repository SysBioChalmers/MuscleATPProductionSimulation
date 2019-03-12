function model = setupSimulation(model, param)
    model = addInternalConstraints(model, param);
    model = addTransportConstraints(model,   's', {'O2'}, -param.vO2max, 0);
    model = addTransportConstraints(model, 'sm3', {'palmitate', 'L-lactate'}, [-param.peripheralFA -param.peripheralLactateCapacity],  [0 0]);    
    model = addTransportConstraints(model, 'sm1', {'palmitate','L-lactate'}, [-1000, -1000],  [0 0]);   
    model = addTransportConstraints(model, 'sm2', {'palmitate','L-lactate'}, [-1000 0], [0 1000]);
    model = addMaintenanceConstraints(model, 1000, param.maintainance, param.internalWork);
    
    model = addExchangeMedum(model);

end

