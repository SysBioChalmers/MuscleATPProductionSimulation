function model = addInternalConstraints(model, param)    

    model = configureSMatrix(model, -1, 'ATPworkOut_m1', 'ATPwork[s]');
    model = configureSMatrix(model, -1, 'ATPworkOut_m2', 'ATPwork[s]');   
    model = configureSMatrix(model, 0, 'ATPworkOut_m1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 0, 'ATPworkOut_m2', 'ATPwork[sm2]');
   
    model = configureSMatrix(model, 0, 'JoinMuscleATP1', 'ATPwork[sm1]');
    model = configureSMatrix(model, 0, 'JoinMuscleATP1', 'ATPwork[sm2]'); 
    model = configureSMatrix(model, 1-param.type2target, 'JoinMuscleATP1', 'ATPwork[cm1]');
    model = configureSMatrix(model, param.type2target, 'JoinMuscleATP1', 'ATPwork[cm2]');     
    
    model.ub(findIndex(model.rxns, 'ATPworkOut_m1')) = param.type1tresh;
    model.ub(findIndex(model.rxns, 'ATPworkOut_m2')) = 0;
    
    scaledM1 = param.dwMuscle  * param.scalingFactor * param.muscleRatio;
    scaledM2 = param.dwMuscle  * param.scalingFactor * (1-param.muscleRatio) * param.m2Efficency;
    
    %Specific activity Complex I
    model.ub(findIndex(model.rxns,'HMR_6921_m1')) = param.HMR_6921 * scaledM1;
    model.ub(findIndex(model.rxns,'HMR_6921_m2')) = param.HMR_6921 * scaledM2;

    %Specific activity Complex IV
    model.ub(findIndex(model.rxns,'HMR_6914_m1')) = param.HMR_6914 * scaledM1;
    model.ub(findIndex(model.rxns,'HMR_6914_m2')) = param.HMR_6914 * scaledM2;

    %Specific activity ETF
    model.ub(findIndex(model.rxns,'HMR_6911_m1')) = param.HMR_6911 * scaledM1;
    model.ub(findIndex(model.rxns,'HMR_6911_m2')) = param.HMR_6911 * scaledM2;

    %Lactate accumulation flux for non steady state
    model.ub(findIndex(model.rxns, 'lactateBuffering')) = param.lactateBuffering;
    
end

