function model = addInternalConstraints(model, dwMuscle, muscleRatio, vO2perDryweight, m2Efficency, complex1Ratio, FAratio, type1tresh, type2target)    

    model = configureSMatrix(model, -1, 'ATPworkOut_m1', 'ATPwork[s]');
    model = configureSMatrix(model, -1, 'ATPworkOut_m2', 'ATPwork[s]');   
    model = configureSMatrix(model, 0, 'ATPworkOut_m1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 0, 'ATPworkOut_m2', 'ATPwork[sm2]');
   
    model = configureSMatrix(model, 0, 'JoinMuscleATP1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 0, 'JoinMuscleATP1', 'ATPwork[sm2]'); 
    model = configureSMatrix(model, 1-type2target, 'JoinMuscleATP1', 'ATPwork[cm1]');
    model = configureSMatrix(model, type2target, 'JoinMuscleATP1', 'ATPwork[cm2]');     
    
    model.ub(findIndex(model.rxns, 'ATPworkOut_m1')) = type1tresh;
    model.ub(findIndex(model.rxns, 'ATPworkOut_m2')) = 0;
    
    scaledC1 = dwMuscle*vO2perDryweight*muscleRatio;
    scaledC2 = dwMuscle*vO2perDryweight*(1-muscleRatio) * m2Efficency;
    

    %the vO2 on complex I substrate vs (I + II)
    %the stochiometric relation to complex IV
    cmpxIfactor = 5/3 * complex1Ratio; 
    
    fattyAcidFactor = FAratio * 1/23; %stochiometery for fatty acids
    
    model.ub(findIndex(model.rxns, 'HMR_6914_m1')) = scaledC1; %Oxygen in m1
    model.ub(findIndex(model.rxns, 'HMR_6914_m2')) = scaledC2; %Oxygen in m2
   
    %complex I
    model.ub(findIndex(model.rxns, 'HMR_6921_m1')) = cmpxIfactor * scaledC1; %Cplx IV m2x
    model.ub(findIndex(model.rxns, 'HMR_6921_m2')) = cmpxIfactor * scaledC2; %Cplx IV m2x
    
    model.lb(findIndex(model.rxns, 'HMR_0217_m1')) = -fattyAcidFactor * scaledC1; %fatty acid metabolism
    model.lb(findIndex(model.rxns, 'HMR_0217_m2')) = -fattyAcidFactor * scaledC2; %fatty acid metabolism  
   
    
end

