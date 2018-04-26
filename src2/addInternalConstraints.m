function model = addInternalConstraints(model, mScale, m1Ratio, c1, c2)
    model = configureSMatrix(model, m1Ratio, 'JoinMuscleATP1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 1-m1Ratio, 'JoinMuscleATP1', 'ATPwork[sm2]');
    
    scaledC1 = mScale*c1*m1Ratio;
    scaledC2 = mScale*c2*(1-m1Ratio);

    model.ub(findIndex(model.rxns, 'HMR_6914_m1')) = scaledC1; %Cplx IV m1
    model.ub(findIndex(model.rxns, 'HMR_6914_m2')) = scaledC2; %Cplx IV m2
    model.ub(findIndex(model.rxns, 'HMR_6921_m1')) = 1/5*5/3*scaledC1; %Cplx IV m1
    model.ub(findIndex(model.rxns, 'HMR_6921_m2')) = 1/5*5/3*scaledC2; %Cplx IV m2
end

