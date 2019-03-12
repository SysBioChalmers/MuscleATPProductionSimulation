function appendCplx1bypass(model, ATPrate, fullSolution)
color1 = [93 93 93]/256;

%Collect simulationdata:
mTransp = getTransport(model, {'O2','CO2'}, 'sb', 's');
mTransp1 = getTransport(model, {'L-lactate'}, 'sb', 'sm2');

O2 = -fullSolution(:,mTransp(1));
CO2 = fullSolution(:,mTransp(2));
lactateFlux = fullSolution(:,mTransp1(1));

O2mod = molToMl(O2);
CO2mod = molToMl(CO2);
RQ = CO2mod./O2mod;
Slactate = convertVlacToConc(lactateFlux, 1.4, 8);
Wmod = molToW(1000*ATPrate);

subplot(2,2,1)
hold all
plot(Wmod, O2mod, 'linewidth', 2, 'color', color1);

subplot(2,2,2)
hold all
plot(Wmod, CO2mod, 'k', 'linewidth', 2, 'color', color1);

subplot(2,2,3)
hold all
plot(Wmod, RQ, 'linewidth', 2, 'color', color1);

subplot(2,2,4)
hold all
plot(Wmod, Slactate, 'linewidth', 2, 'color', color1);
end
