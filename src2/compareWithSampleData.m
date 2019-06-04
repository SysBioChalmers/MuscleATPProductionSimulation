function compareWithSampleData(model, dataFolder, ATPrate, fullSolution)
color3 = [204 204 204]/256;
color2 = [215 86 40]/256;
color1 = [93 155 211]/256;

%Collect simulationdata:
mTransp = getTransport(model, {'O2','CO2'}, 'sb', 's');
O2 = -fullSolution(:,mTransp(1));
CO2 = fullSolution(:,mTransp(2));

O2mod = molToMl(O2);
CO2mod = molToMl(CO2);
RQ = CO2mod./O2mod;
RQtresh = find(RQ>0.98);
if not(isempty(RQtresh))
    RQtresh = RQtresh(1);
else
    RQtresh = length(ATPrate);
end

Wmod = molToW(1000*ATPrate);


complexIbypassM1 = getComplexIbypass(model, fullSolution, 'HMR_6921_m1');
complexIbypassM2 = getComplexIbypass(model, fullSolution, 'HMR_6921_m2');

%Collect sample data
data = loadSampleData(dataFolder);
plotExpData(data, color1)
% data = importdata(['sampleData/' dataFolder '/data2.txt']);
% plotExpData(data, color1)
%  data = importdata(['sampleData/' dataFolder '/data3.txt']);
%  plotExpData(data, color1)

subplot(2,2,1)
hold all
plot(Wmod, O2mod, 'linewidth', 2, 'color', color2);
ylim([0 inf])
%set(gca,'XTick',[])
%ylabel('ml/min')
ylabel('vO2 [ml/min]')

% ivalues = round(length(Wmod)*[6/20 14/20]);
% intervals = round(4/20*length(Wmod));
% 
% for i = ivalues
%     deltaEff = (O2mod(i+intervals)-O2mod(i))/(Wmod(i+intervals)-Wmod(i));
%     area([Wmod(i) Wmod(i+intervals)],  [O2mod(i) O2mod(i+intervals)], 'FaceColor', color3, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
%     text(Wmod(i), 300, sprintf('\\Delta%2.1f',deltaEff), 'Interpreter', 'tex');
% end
%plot(Wmod(RQtresh) * [1 1], [0 4500],'k--')
%plot(Wmod(complexIbypassM2) * [1 1], [0 4500],'k--')
area([Wmod(RQtresh) Wmod(complexIbypassM2)],  [4500 4500], 'FaceColor', color3, 'FaceAlpha', 0.5, 'EdgeColor', 'none')


subplot(2,2,2)
hold all
plot(Wmod, CO2mod, 'k', 'linewidth', 2, 'color', color2);
ylabel('vCO2 [ml/min]')
ylim([0 inf])
%plot(Wmod(RQtresh) * [1 1], [0 5000],'k--')
%plot(Wmod(complexIbypassM2) * [1 1], [0 5000],'k--')
area([Wmod(RQtresh) Wmod(complexIbypassM2)],  [5000 5000], 'FaceColor', color3, 'FaceAlpha', 0.5, 'EdgeColor', 'none')

subplot(2,2,3)
hold all
plot(Wmod, RQ, 'linewidth', 2, 'color', color2);
ylabel('RER (RQ)')
xlabel('W')
ylim([0.7 inf])

plot(Wmod(RQtresh) * [0 1], [1 1],'k--')
%plot(Wmod(RQtresh) * [1 1], [0 1.11],'k--')
%plot(Wmod(complexIbypassM2) * [1 1], [0 1.11],'k--')
area([Wmod(RQtresh) Wmod(complexIbypassM2)],  [1.11 1.11], 'FaceColor', color3, 'FaceAlpha', 0.5, 'EdgeColor', 'none')


subplot(2,2,4)
hold all
mTransp1 = getTransport(model, {'L-lactate'}, 'sb', 'sm2');
%mTransp2 = getTransport(model, {'L-lactate'}, 'sb', 'sm3');

lactateFlux = fullSolution(:,mTransp1(1));
if isfile(['sampleData/' dataFolder '/la.txt'])
    lactate = load(['sampleData/' dataFolder '/la.txt']);
else
    lactate = zeros(1,2);
end
%plot(lactate(:,1), lactate(:,2), 'o--');

Slactate = convertVlacToConc(lactateFlux, 1.4, 8);

ylim([0 inf])
plot(Wmod, Slactate, 'linewidth', 2, 'color', color2);

%plot(Wmod(RQtresh) * [1 1], [0 15],'k--')
%plot(Wmod(complexIbypassM2) * [1 1], [0 15],'k--')
area([Wmod(RQtresh) Wmod(complexIbypassM2)],  [15 15], 'FaceColor', color3, 'FaceAlpha', 0.5, 'EdgeColor', 'none')


scatter(lactate(:,1),lactate(:,2), 15, 'MarkerFaceColor', color1,'MarkerEdgeColor', color1, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
xlabel('W')
ylabel('lactate [mM]')
%ylabel('lactate mM')

% figure()
% hold on
% plotScatter(Wdata, Wdata./VO2, color1)  
% plot(Wmod, Wmod./O2mod, 'linewidth', 2, 'color', color2);
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
   fill(xvals, yvals, color, 'edgecolor','none', 'FaceAlpha', 0.3)
end

function plotScatter(X, Y, color)  
   scatter(X,Y, 10, 'MarkerFaceColor', color,'MarkerEdgeColor', color, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
end

% function VE = ventilationModel(gasExchange)
%     SGas = revMM(gasExchange', 5000, 8000);
%     VE = MM(SGas, 180, 17000);
% end


function plotExpData(data, color)
    Wdata = data(:,1);
    VO2 = data(:,2);
    VCO2 = data(:,3);
    RQest = VCO2./VO2;
    wMax = max(Wdata)
    
    k = 8;
    tickValues = 0:50:wMax;

    subplot(2,2,1)
    hold all
    X = movmean(Wdata,k);
    Y = movmean(VO2,k);
    E = movstd(VO2,k);
    %plotErrorArea(X, Y, E, color)
    plotScatter(Wdata, VO2, color)  
    xticks(tickValues)
    xlim([0 wMax])
    
    subplot(2,2,2)
    hold all
    Y = movmean(VCO2,k);
    E = movstd(VCO2,k);
    %plotErrorArea(X, Y, E, color)  
    plotScatter(Wdata, VCO2, color)
    xticks(tickValues)
    xlim([0 wMax])

    subplot(2,2,3)
    hold all
    Y = movmean(RQest,k);
    E = movstd(RQest,k);
    %plotErrorArea(X, Y, E, color)  
    plotScatter(Wdata, RQest, color)
    xticks(tickValues)
    xlim([0 wMax])

    subplot(2,2,4)
    hold all
    xticks(tickValues)
    xlim([0 wMax])  
end

