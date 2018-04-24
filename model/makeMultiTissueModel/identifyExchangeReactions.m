function [allReactions, metaboliteNR] = identifyExchangeReactions(modelS)
    elementExists = abs(sign(modelS));
    allReactions = find(sum(elementExists)==1);
    metaboliteNR = zeros(length(allReactions), 1);
    for j = 1:length(allReactions)
        metaboliteNR(j) = find(elementExists(:,allReactions(j)));
    end
end

