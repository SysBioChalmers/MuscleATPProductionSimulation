function smoothData = moveAvg(dataX, dataY, lag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

bothData = [dataX, dataY];
deletedData = or(isnan(bothData(:,1)), isnan(bothData(:,2)));
bothData(deletedData,:) = [];

[crap, orderData] = sort(bothData(:,1));
bothData = bothData(orderData,:);

smoothData = tsmovavg(bothData,'s',lag,1);
deletedData = or(isnan(smoothData(:,1)), isnan(smoothData(:,2)));
smoothData(deletedData,:) = [];

end

