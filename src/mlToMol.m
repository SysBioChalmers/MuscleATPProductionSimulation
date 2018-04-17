function result = mlToMol(data)
    %Converts ml/min to mMol/h
    %1 ?mol O2 = .022391 ml
    result = data/(0.022391 * 10^6);
    result = result * 60;


end

