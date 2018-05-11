function Slactate = convertVlacToConc(lactateFlux, baseLine, vmax)
%    vmax = 9;
    km = 10.73;
%    baseLine = 1.3; 
    basalFlux = MM(baseLine, vmax, km);
    lactateFlux = lactateFlux + basalFlux;
    Slactate = revMM(lactateFlux, vmax, km);
    %https://www.sigmaaldrich.com/catalog/product/roche/lldhro?lang=en&region=SE
    
end

function v = MM(S, vmax,km)
    v = vmax*S./(km+S);
end

function S = revMM(v,vmax,km)
    S = (km*v)./(vmax-v);
end
