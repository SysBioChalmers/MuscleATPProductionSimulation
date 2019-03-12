function param = fitModel(model, param, fitList, subjectName)
    N = 30;
    range = 0.6;

    [allData, vO2max] = loadData(subjectName);
    settings = setSimulationSettings(N, false);
     
    param.vO2max = mlToMol(vO2max);
    
    dataX = allData(:,1);
    dataY = allData(:,2:3);
    
    %Reactions to plot
    readouts = getTransport(model, {'O2', 'CO2'}, 'sb', 's');

    %Set intial point
    x0 = zeros(length(fitList),1);
    lb = zeros(length(fitList),1);
    ub = zeros(length(fitList),1);
    for i = 1:length(fitList)
        x0(i) = param.(fitList{i});
        lb(i) = param.(fitList{i})*(1-range);
        ub(i) = param.(fitList{i})*(1+range);        
    end    
    
    opts = optimset('Display','off');
    [p,resnormMin,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@(x,y) fitFunction(x,y, model, settings, param, fitList, readouts),x0,dataX,dataY,lb,ub,opts);

    
    for i = 1:length(fitList)
        param.(fitList{i}) = p(i);
    end    
    
    %conf = nlparci(x,residual,'jacobian',jacobian);
end


function y = fitFunction(p, x, model, settings, param, fitList, readouts)
    disp('iter')
    
    for i = 1:length(fitList)
        param.(fitList{i}) = p(i);
    end
        
    model = setupSimulation(model, param);
    [ATPrate, fullSolution] = runFullModel(model, settings);

    if size(fullSolution,1) == 1
        size(ATPrate)
        N = settings.timeSteps;
        ATPrate = linspace(0,100,N);
        fullSolution = rand(N,length(model.rxns));
    end
    
    Wmod = molToW(1000*ATPrate)';    
    
    vO2mol = -fullSolution(:, readouts(1));
    vO2 = interp1q(Wmod,molToMl(vO2mol),x);
    vCO2mol = fullSolution(:, readouts(2));
    vCO2 = interp1q(Wmod,molToMl(vCO2mol),x);
    
    y = [vO2 vCO2];
end


function [allData, vO2max] = loadData(dataFolder)
    if strcmp(dataFolder, 'subject1')
        %temp for subject 1
        data = importdata(['sampleData/' dataFolder '/data1.txt']);
        dataHigh = data.data;
        dataHigh(dataHigh(:,1)>300,:) = [];
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
