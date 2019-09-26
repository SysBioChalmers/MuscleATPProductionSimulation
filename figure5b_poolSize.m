

symbolicTime = linspace(1, 20);

pool = 1;

flux = pool./symbolicTime;
symbolicTime = [0 symbolicTime];
flux = [flux(1) flux];

hold all
plot(symbolicTime, flux)
area(symbolicTime, flux)


ylim([0 0.5])
xlim([0 20])
