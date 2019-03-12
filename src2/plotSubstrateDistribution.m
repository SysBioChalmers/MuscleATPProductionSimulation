function plotSubstrateDistribution(model, ATPrate, fullSolution)
roundedWatt = round(1000*molToW(ATPrate));

metList = {'glycogen', 'L-lactate', 'palmitate'};
carbonContent = [6 3 16];

otherComp={'sm1', 'sm2'};

wMax = max(roundedWatt);
wMax = 372;

tickValues = 0:50:wMax;

%colors = get(gca,'colororder');
for i = 1:2
    subplot(1,2,i)
    transportReactions = getTransport(model, metList, 'sb', otherComp{i});
    fluxes = -fullSolution(:, transportReactions);
    normalizedFlux = fluxes.*repmat(carbonContent, length(ATPrate),1);
    normalizedFlux(normalizedFlux<0) = 0;
    totalCarbon = sum(normalizedFlux,2);
    size(totalCarbon)
    Yvalues = normalizedFlux./repmat(totalCarbon, 1,length(metList));
    area(roundedWatt, Yvalues, 'EdgeColor', 'none')
    xticks(tickValues)
    xlim([0 wMax]);
    ylim([0 1])
    xlabel('W')
end
legend(metList)
legend boxoff
end

