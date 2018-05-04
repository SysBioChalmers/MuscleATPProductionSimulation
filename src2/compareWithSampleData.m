function compareWithSampleData(model, dataFolder, ATPrate, fullSolution, lag)
%Expected baseline ~ 245 ml O2/min

factorForInternalWork = 0;

%Collect simulationdata:
mTransp = getTransport(model, {'O2','CO2'}, 'sb', 's');
O2 = -fullSolution(:,mTransp(1));
CO2 = fullSolution(:,mTransp(2));
O2 = molToMl(O2);
CO2 = molToMl(CO2);
W = molToW(1000*ATPrate);

CO2mod = tsmovavg(CO2','s',lag,2);
O2mod = tsmovavg(O2','s',lag,2);
Wmod = tsmovavg(W,'s',lag,2);
% CO2mod = CO2;
% O2mod = O2;
% Wmod = W; 

ventMod = max(O2mod, CO2mod);
ventMod = ventilationModel(max(O2mod, CO2mod));%43 Ml per breath

%Collect sample data
% VO2 = load(['sampleData/' data '/VO2.txt']);
% VCO2 = load(['sampleData/' data '/VCO2.txt']);
% RQest = load(['sampleData/' data '/RQ.txt']);
% ventData = load(['sampleData/' data '/VE.txt']);
data = importdata(['sampleData/' dataFolder '/data.txt']);
data = data.data;
Wdata = data(:,1);
VO2 = data(:,3);
VCO2 = data(:,4);
ventData = data(:,5);
RQest = VCO2./VO2;


%Fit polymodel:
%Hmodel = polyfit(HRData(:,1), HRData(:,2),2);
%HRestimate = Hmodel(1)*(O2mod).^2 + Hmodel(2)*(O2mod) + Hmodel(3);

subplot(2,2,1)
hold all
k = 8;
X = movmean(Wdata,k);
Y = movmean(VO2,k);
E = movstd(VO2,k);
color1 = [250 50 80]/256;
color2 = [80 50 250]/256;
plotErrorArea(X, Y, E, color2)
scatter(Wdata,VO2, 10, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.8)

plot(Wmod, O2mod, 'k', 'linewidth', 2);

ylabel('ml/min')
h = findobj(gca,'Type','line');

ylim([0 inf])
set(gca,'XTick',[])
ylabel('ml/min')

subplot(2,2,2)
hold all
Y = movmean(VCO2,k);
E = movstd(VCO2,k);

plotErrorArea(X, Y, E, color2)  
scatter(Wdata,VCO2, 10, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.8)

plot(Wmod, CO2mod, 'k', 'linewidth', 2);
%scatter(VCO2(:,1), VCO2(:,2), 'bo');
%scatter(VO2(:,1), VO2(:,2), 'ro');


ylabel('ml/min')
set(gca,'XTick',[])
h = findobj(gca,'Type','line');
ylim([0 inf])

subplot(2,2,3)
hold all
Y = movmean(RQest,k);
E = movstd(RQest,k);
plotErrorArea(X, Y, E, color2)  
scatter(Wdata,RQest, 10, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.8)

plot(Wmod, CO2mod./O2mod, 'k', 'linewidth', 2);

%scatter(RQest(:,1), RQest(:,2), 'bo');
ylabel('RQ')
%set(gca,'XTick',[])
% subplot(2,2,3)
% hold all
% Y = movmean(ventData,k);
% E = movstd(ventData,k);
% plotErrorArea(X, Y, E, color2)  
% 
% 
% plot(Wmod, ventMod, 'linewidth', 2, 'color', color2);
% scatter(Wdata, ventData, 10, 'MarkerFaceColor', color2, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
% 
% xlabel('W')
% ylabel('VE')
% 
% subplot(2,2,4)
% hold all
% 
% delta = diff([Wmod' O2mod']);
% efficency = delta(:,1)./delta(:,2);
% %efficency = Wmod./O2mod;
% plot(Wmod(2:end), efficency, 'linewidth', 3);
% 
% %timePoints = linspace(min(VO2(:,1)), max(VO2(:,1)), 30);
% %linReg = polyfit(VO2(:,1),VO2(:,2),2);
% %plot(timePoints, timePoints*linReg(1)*2 + linReg(2), 'k-', 'linewidth', 2)
% 
% 
% efficency2 = diff(smoothO2);
% efficency2 = efficency2(:,1)./efficency2(:,2);
% 
% plot(smoothO2(2:end,1), efficency2, 'bo')
% trend = tsmovavg([smoothO2(2:end,1) efficency2],'s',2,1);
% plot(trend(:,1), trend(:,2), 'k-', 'linewidth', 2)
% maxVal = max(max(efficency), max(efficency2))*1.1;
% ylim([0 maxVal])
% 
% xlabel('W')
% %ylabel('dW/dO2')
% ylabel('D Efficency')

%legend('Model', 'Data (moving avg.)', 'trend (moving avg.)', 'location', 'se')


%plot(O2mod, HRestimate, 'linewidth', 3);
%scatter(HRData(:,1), HRData(:,2), 'bx');
%xlabel('vO2')
%ylabel('Heart rate')


subplot(2,2,4)
hold all
mTransp = getTransport(model, {'L-lactate'}, 'sb', 'sm2');
lactateFlux = fullSolution(:,mTransp(1));
lactateMod = tsmovavg(lactateFlux,'s',lag,1);
lactate = load(['sampleData/' dataFolder '/la.txt']);
%plot(lactate(:,1), lactate(:,2), 'o--');
scatter(lactate(:,1), lactate(:,2), 10, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.8)


baseLine = mean(lactate(1:4,2));

vmax = 9;
km = 10.73;
%https://www.sigmaaldrich.com/catalog/product/roche/lldhro?lang=en&region=SE

basalFlux = MM(baseLine, vmax, km);

lactateFlux = lactateMod + basalFlux;
Slactate = revMM(lactateFlux, vmax, km);
%lactateFlux = lactateFlux 
%lactateConc = lactateConc + baseLine;

ylim([0 inf])
plot(W, Slactate, 'k-', 'linewidth', 2);

xlabel('W')
ylabel('lactate mM')


end

function v = MM(S, vmax,km)
    v = vmax*S./(km+S);
end

function S = revMM(v,vmax,km)
    S = (km*v)./(vmax-v);
end

function plotErrorArea(X, Y, E, color)  

   lb = Y-E;
   ub = Y+E;
   
   xvals = [X; flipud(X)];
   yvals = [lb; flipud(ub)];  
   fill(xvals, yvals, color, 'edgecolor','none', 'FaceAlpha',0.4)
   
end

function VE = ventilationModel(gasExchange)
    SGas = revMM(gasExchange', 5000, 8000);
    VE = MM(SGas, 180, 17000);
end


