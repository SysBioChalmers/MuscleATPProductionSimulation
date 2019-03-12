function [delta, se, p, vO2max, intercept, RQout] = fitSubject(subjectName, O2corr)
[dataLow, dataHigh, vO2max] = loadData(subjectName, O2corr);
allData = [dataLow;dataHigh];

RQ = allData(:,3)./allData(:,2);

%intercept = estimateIntercept(dataLow, dataHigh, midLower, midUpper);
intercept = mixedModelIntercept(allData);
midpoint = intercept(1,1);


W = allData(:,1)-midpoint;
vO2 = allData(:,2);
LH = W>0;
RQout = [median(RQ(not(LH))) median(RQ(LH))];

data = table(W, vO2, LH);
mdlLow = fitlm(data, 'vO2~W:(1+LH)');
p = mdlLow.Coefficients.pValue(3);

LH = (LH==0);
data = table(W, vO2, LH);
mdlHigh = fitlm(data, 'vO2~W:(1+LH)');
delta = [mdlLow.Coefficients.Estimate(2);mdlHigh.Coefficients.Estimate(2)];
se = [mdlLow.Coefficients.SE(2);mdlHigh.Coefficients.SE(2)];

%intercept(2,1)/vO2max

ylim([0 5000])
plotRawdata(dataLow, dataHigh, subjectName);
plotIntercept(intercept);
plotFits(mdlHigh, data, midpoint);
printText(delta, se, intercept,  p);
plot([0 300], vO2max * [1 1], '-'); 
xlabel('W')
ylabel('vO2') 
if O2corr
    ylabel('RQ adjusted vO2')
end
xlim([0 300])
end



function [dataLow, dataHigh, vO2max] = loadData(dataFolder, O2corr)
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


    if O2corr == true
        dataLow(:,2) = normalizeO2(dataLow(:,2),dataLow(:,3));
        dataHigh(:,2) = normalizeO2(dataHigh(:,2),dataHigh(:,3));
    end

    data = importdata(['sampleData/' dataFolder '/data1.txt']);
    vO2max = max(data.data(:,2));
end

function intercept = estimateIntercept(dataLow, dataHigh, midLower, midUpper)
    resolution = 2000;

    %separate data into low and high
    dataLow = dataLow(dataLow(:,1)<=midLower,:);
    dataHigh = dataHigh(dataHigh(:,1)>=midUpper,:);

    mdl1 = fitlm(dataLow(:,1), dataLow(:,2));
    mdl2 = fitlm(dataHigh(:,1), dataHigh(:,2));

    mdl1k = mdl1.Coefficients.Estimate;
    mdl2k = mdl2.Coefficients.Estimate;

    interceptX = (mdl2k(1)-mdl1k(1))/(mdl1k(2)-mdl2k(2));
    interceptY = mdl1k(1) + interceptX*mdl1k(2);

    xvalues=linspace(0, max(dataHigh(:,1)),resolution);

    [ypred1,yci1] = predict(mdl1, xvalues');
    [ypred2,yci2] = predict(mdl2, xvalues');

    overlap = sum((yci1-yci2)>0,2)==1;
    
    overlapX = [interceptX xvalues(overlap)];
    overlapY = [interceptY yci1(overlap,1)' yci2(overlap,2)'];
    
    intercept = [
        interceptX min(overlapX) max(overlapX)
        interceptY min(overlapY) max(overlapY)
        ];
end

function y = fitFunction(p,x)
    y = zeros(size(x,1),1);
    
    crossOver = p(4);
    delta1 = p(2);
    delta2 = p(3);
    m1 = p(1);
    m2 = m1 + delta1*crossOver;
        
    lowX = x<=crossOver;
    highX = x>crossOver;
    
    y(lowX) = m1 + delta1*x(lowX);
    y(highX) = m2 + delta2*(x(highX)-crossOver);
end


function intercept = mixedModelIntercept(allData)
    dataX = allData(:,1);
    dataY = allData(:,2);
    randomRestarts = 100;
    
    x0(1) = 800;
    x0(2) = 9;
    x0(3) = 12;
    x0(4) = 150;
    opts = optimset('Display','off');
    resnormMin = inf;
    for i = 1:randomRestarts
        x0(4) = x0(4) + (rand-0.5)*200;
        [curX,curResnorm,curResidual,exitflag,output,lambda,curJacobian] = lsqcurvefit(@(x,y) fitFunction(x,y),x0,dataX,dataY,[],[],opts);
        if curResnorm<resnormMin
            x = curX;
            resnormMin = curResnorm;
            residual = curResidual;
            jacobian = curJacobian;
        end
    end
    
    conf = nlparci(x,residual,'jacobian',jacobian);

    interceptX = [x(4) conf(4,1) conf(4,2)];
    interceptY = fitFunction(x,interceptX);

    intercept = [
        interceptX
        interceptY
        ];
end

function plotRawdata(dataLow, dataHigh, subjectName)
    hold all
    color1 = [93 155 211]/256;
    color2 = [215 86 40]/256;

%     [ypred1,yci1] = predict(mdlLow, [xvals' xvals'>0],'alpha',0.01);
%     [ypred2,yci2] = predict(mdlLow, [xvals' xvals'>inf]);

    %plot rawdata:
    scatter(dataLow(:,1), dataLow(:,2), 15, 'filled', 'markerfacecolor', color1)
    scatter(dataHigh(:,1), dataHigh(:,2), 15, 'filled', 'markerfacecolor', color2)
    legend({'session-low', 'session-high'}, 'location', 'nw')
    legend boxoff
    title(subjectName)
end

function plotIntercept(intercept)
    %plot intercept
%    ymax = ylim;
%     plot(midLower * [1 1], [0 ymax(2)], 'k--','HandleVisibility','off');
%     plot(midUpper * [1 1], [0 ymax(2)], 'k--','HandleVisibility','off');
    xvals = intercept(1,:);
    yvals = intercept(2,:);
    
    fillX = [0 xvals(2) xvals(2) xvals(3)];
    fillY = [[yvals(2);yvals(2);0;0] [yvals(3);yvals(3);yvals(3);yvals(3)]];
    fillArea(fillX, fillY, [0.6 0.6 0.6])
    plot(xvals(1) * [0 1 1], yvals(1) * [1 1 0], '-', 'color', [0.5 0.5 0.5], 'HandleVisibility','off')
end

function plotFits(mdlHigh, data, midpoint)
    xvalues=linspace(min(data.W),max(data.W),100);
    xIn = [xvalues' xvalues'<0];    
    [ypred1,yci1] = predict(mdlHigh, xIn);

    xIn = [xvalues' ones(length(xvalues),1)];    
    [ypred2,yci2] = predict(mdlHigh, xIn);    
    xIn = [xvalues' zeros(length(xvalues),1)];    
    [ypred3,yci3] = predict(mdlHigh, xIn);
    
    xvalues = xvalues + midpoint;
    
    fillArea(xvalues, yci1, [0.4 0.4 0.4])
    plot(xvalues, ypred2, 'k--', 'linewidth', 2, 'HandleVisibility','off');
    plot(xvalues, ypred3, 'k--', 'linewidth', 2, 'HandleVisibility','off');
    plot(xvalues, ypred1, 'k-', 'linewidth', 2, 'HandleVisibility','off');
    
%     fillArea(xvalues, yci2, [0.4 0.4 0.4])
%     plot(xvalues, ypred2, 'k-', 'linewidth', 2,'HandleVisibility','off');    
end


function printText(delta, se, intercept,  p) 
    %Make text for slopes:
    text(intercept(1,3), 1000, sprintf('\\Delta1=%2.1f\\pm%2.1f', delta(1), se(1)))
    text(intercept(1,3), 1500, sprintf('\\Delta2=%2.1f\\pm%2.1f', delta(2), se(2)))
    text(intercept(1,3), 500, sprintf('p=%2.2e', p))
end