function model = addInternalConstraints(model, m1, m2)
    model.ub(findIndex(model.rxns, 'HMR_6916_m1')) = m1; %ATP synthase m1
    model.ub(findIndex(model.rxns, 'HMR_6916_m2')) = m2; %ATP synthase m1
    %Test Knock out muscle komplex 1
    %model.ub(findIndex(model.rxns, 'HMR_6921_m1')) = 0;
    %model.ub(findIndex(model.rxns, 'HMR_6921_m2')) = 0;
end

