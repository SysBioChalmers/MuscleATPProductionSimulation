function fullSolution = crossAndScale(fullSolution, ATPrate, low, high, intensity)
    mixLow = find(ATPrate>low);
    mixLow = mixLow(1);
    mixHigh = find(ATPrate>high);
    mixHigh = mixHigh(1);    
    mixpoints = mixLow:mixHigh;
    mixAmount = linspace(0,1, length(mixpoints));
    mixMatrix = repmat(mixAmount',1, size(fullSolution,2));
    A = repmat(fullSolution(mixHigh,:), length(mixpoints), 1).*mixMatrix;
    B = repmat(fullSolution(mixLow,:), length(mixpoints), 1).*(1-mixMatrix);
    fullSolution(mixpoints,:) = intensity * (A+B) + (1-intensity) * fullSolution(mixpoints,:);
end

