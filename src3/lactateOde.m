function [dydt] = lactateOde(t, y, volume, basalFlux, vmax, muscleflux)
    S = y(1)/volume; %mM lactate
    v1 = 1000*basalFlux; %mM/h
    v2 = 1000*muscleflux;     %mM/h
    v3 = 1000*mmLactate(S, vmax); %mM/h
    
    dydt = (v1 + v2 - v3);
end

