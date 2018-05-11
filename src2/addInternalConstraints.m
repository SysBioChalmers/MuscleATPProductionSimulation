function model = addInternalConstraints(model, dwMuscle, m1Ratio, vO2perDryweight, m2Efficency, complex1Ratio, FAratio)
    model = configureSMatrix(model, m1Ratio, 'JoinMuscleATP1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 1-m1Ratio, 'JoinMuscleATP1', 'ATPwork[sm2]');
    
    scaledC1 = dwMuscle*vO2perDryweight*m1Ratio;
    scaledC2 = dwMuscle*m2Efficency*vO2perDryweight*(1-m1Ratio);

    %the vO2 on complex I substrate vs (I + II)
    %the stochiometric relation to complex IV
    cmpxIfactor = 5/3 * complex1Ratio; 
    
    
    fattyAcidFactor = FAratio * 1/23; %stochiometery for fatty acids
    
    model.ub(findIndex(model.rxns, 'HMR_6914_m1')) = scaledC1; %Oxygen in m1
    model.ub(findIndex(model.rxns, 'HMR_6914_m2')) = scaledC2; %Oxygen in m2
    model.ub(findIndex(model.rxns, 'HMR_6921_m1')) = cmpxIfactor * scaledC1; %Cplx IV m1
    model.ub(findIndex(model.rxns, 'HMR_6921_m2')) = cmpxIfactor * scaledC2; %Cplx IV m2
    
    model.lb(findIndex(model.rxns, 'HMR_0217_m1')) = -fattyAcidFactor * scaledC1; %fatty acid metabolism
    model.lb(findIndex(model.rxns, 'HMR_0217_m2')) = -fattyAcidFactor * scaledC2; %fatty acid metabolism  
    
    
end

