function model = addInternalConstraints(model, dwMuscle, m1Ratio, c1, c2)
    model = configureSMatrix(model, m1Ratio, 'JoinMuscleATP1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 1-m1Ratio, 'JoinMuscleATP1', 'ATPwork[sm2]');
    
    scaledC1 = dwMuscle*c1*m1Ratio;
    scaledC2 = dwMuscle*c2*(1-m1Ratio);

    cmpxIfactor = 33.45/76.41; %the vO2 on complex I substrate vs (I + II)
    cmpxIfactor = 5/3 * cmpxIfactor; %the stochiometric relation to complex IV
    
    fattyAcidFactor = 16.55/76.41;
    fattyAcidFactor = fattyAcidFactor * 1/23; %stochiometery for fatty acids
    fattyAcidFactor = fattyAcidFactor * 0.6; %fitting parameter
    
    model.ub(findIndex(model.rxns, 'HMR_6914_m1')) = scaledC1; %Oxygen in m1
    model.ub(findIndex(model.rxns, 'HMR_6914_m2')) = scaledC2; %Oxygen in m2
    model.ub(findIndex(model.rxns, 'HMR_6921_m1')) = cmpxIfactor * scaledC1; %Cplx IV m1
    model.ub(findIndex(model.rxns, 'HMR_6921_m2')) = cmpxIfactor * scaledC2; %Cplx IV m2
    
    model.lb(findIndex(model.rxns, 'HMR_0217_m1')) = -fattyAcidFactor * scaledC1; %fatty acid metabolism
    model.lb(findIndex(model.rxns, 'HMR_0217_m2')) = -fattyAcidFactor * scaledC2; %fatty acid metabolism  
    
    
end

