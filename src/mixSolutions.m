function AandB = mixSolutions(A, B, startPoint, stopPoint)
    scaleSize = size(A,1);
    startId = round(startPoint * scaleSize);
    stopId = round(stopPoint * scaleSize);
    mixingPoints = stopId - startId;
    stopLength = scaleSize- stopId;
    mixVector = [zeros(1,startId), linspace(0,1, mixingPoints), ones(1,stopLength)];
    mixMatrix = repmat(mixVector', 1, size(A,2));
    AandB = A.*(1-mixMatrix) + B.*mixMatrix;
end

