function Slactate = convertVlacToConc(lactateFlux, baseConcentration, vmax)
%    vmax = 9;
    km = 10.73;
%    baseLine = 1.3; 
    basalFlux = MM(baseConcentration, vmax, km);
    lactateFlux = lactateFlux + basalFlux;
    Slactate = revMM(lactateFlux, vmax, km);
    Slactate((lactateFlux+basalFlux)>vmax) = inf;
    
    %https://www.sigmaaldrich.com/catalog/product/roche/lldhro?lang=en&region=SE
    
end

function v = MM(S, vmax,km)
    v = vmax*S./(km+S);
end

function S = revMM(v,vmax,km)
    S = (km*v)./(vmax-v);
end
