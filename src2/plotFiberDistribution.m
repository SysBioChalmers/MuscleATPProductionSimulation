function plotFiberDistribution(model, ATPrate, fullSolution, muscleRatio)
    exMap = [67 116 160
             80 137 188
             91 155 213
             151 185 224
             190 209 234]/255;
    hold all
    stochiometry = 0.45;
    normalizedATP = fullSolution(:, findIndex(model.rxns, 'ATPworkOut_m1')) + stochiometry * fullSolution(:, findIndex(model.rxns, 'JoinMuscleATP1'));
    normalizedATP =[normalizedATP./ATPrate' (ATPrate'-normalizedATP)./ATPrate'];
    h = area(ATPrate, normalizedATP, 'FaceColor', exMap(1,:), 'EdgeColor', 'none');
    h(2).FaceColor  = exMap(2,:);
    legend('type1', 'type2')
    plot([min(ATPrate), max(ATPrate)], muscleRatio*[1 1], 'k--', 'linewidth', 2);
    ylim([0 1])
end

