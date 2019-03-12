function [allData, vO2max] = loadSampleData(dataFolder)
    if strcmp(dataFolder, 'subject1')
        %temp for subject 1
        data = importdata(['sampleData/' dataFolder '/data1.txt']);
        dataHigh = data.data;
        dataLow = zeros(0,3);
    else
        data = importdata(['sampleData/' dataFolder '/data3.txt']);
        dataLow = data.data;
        dataLow(dataLow(:,1)<=0,:) = [];
        data = importdata(['sampleData/' dataFolder '/data2.txt']);
        dataHigh = data.data;
    end    
    allData = [dataLow;dataHigh];
    
    data = importdata(['sampleData/' dataFolder '/data1.txt']);
    vO2max = max(data.data(:,2));
end
