function model = addSpecificActivityConstraint(model, cParam, bondParam, time)
    constraintRow = length(model.mets) + 1;

    model.mets{constraintRow} = 'MassConstraint';
    model.metNames{constraintRow} = 'Mass Constraint';
    model.metComps(constraintRow) = 1;
    model.inchis{constraintRow} = '';
    model.metFormulas{constraintRow} = '';
    model.metMiriams{constraintRow} = '';
    model.b = [model.b -model.b];
    
    model.b(constraintRow,2) = bondParam;
    for i = 1:length(model.specificActivity)
        curVal = model.specificActivity(i);
        if curVal ~= -1
            model.S(constraintRow, i) = 1/(time * cParam * curVal);
        end
    end


end

