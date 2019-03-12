function result = mlToMol(data)
    %Converts ml/min to mol/h 
    %22.391 mol/l
    result = data*10^-3/(22.391);
    result = result * 60;
end

