function result = linReg(dataX,dataY,polyFactor)
    if nargin == 2
        polyFactor = 1;
    end
    error = isnan(dataY);
    dataX(error) = [];
    dataY(error) = [];
    xValues = linspace(min(dataX), max(dataX), 100);
    [km, error] = polyfit(dataX, dataY, polyFactor);
    
    yValues = zeros(1,100);
    for i = 1:(polyFactor+1)
        sqr = polyFactor + 1 - i;
        yValues = yValues  + xValues.^sqr * km(i);
    end
    
    result = [xValues' yValues'];
end

