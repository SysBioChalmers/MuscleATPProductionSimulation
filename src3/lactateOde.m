function [dydt] = lactateOde(t, y, volume, basalFlux, muscleflux)
    S = y(1); %mM lactate
    v1 = basalFlux;
    v2 = muscleflux;
    v3 = mmLactate(S);
    
    dydt = (v1 + v2 - v3)/volume;
end
