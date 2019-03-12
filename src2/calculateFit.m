function rmse = calculateFit(referenceData, sampleData)
    punishmentForNan = 1.1;
    interpolData = interp1q(sampleData(:,1), sampleData(:,2), referenceData(:,1));
    fail = isnan(interpolData);
    interpolData(fail) = referenceData(fail,2)*punishmentForNan;
    rmse = sqrt(mean((interpolData-referenceData(:,2)).^2));
    rmse = rmse/mean(referenceData(:,2));
end

