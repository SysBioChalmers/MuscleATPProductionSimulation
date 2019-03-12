hold all
color2 = [215 86 40]/256;
color1 = [93 155 211]/256;
addpath('../src2')
dataFolder = 'subject1';
data = importdata(['../sampleData/' dataFolder '/data1.txt']);
data = data.data;

Wdata = data(:,1);
VO2 = data(:,3);
VCO2 = data(:,4);
EE = O2toenergy(VO2, VCO2);
plot(Wdata, EE,'k.')
xlabel('W')
ylabel('kcal')

normalizedW = Wdata/max(Wdata);

midpoint = 0.45;

intervall1 = and(normalizedW<midpoint, normalizedW>0.1);
intervall2 = and(normalizedW>midpoint, normalizedW<0.9);

xvalues = linspace(0, max(Wdata));

mdl1 = fitlm(Wdata(intervall1), EE(intervall1));
mdl2 = fitlm(Wdata(intervall2), EE(intervall2));

mdl1m = mdl1.Coefficients.Estimate(1);
mdl1k = mdl1.Coefficients.Estimate(2);
mdl2m = mdl2.Coefficients.Estimate(1);
mdl2k = mdl2.Coefficients.Estimate(2);

plot(xvalues, mdl1m + xvalues * mdl1k, 'color', color1, 'linewidth', 2);
plot(xvalues, mdl2m + xvalues * mdl2k, 'color', color2, 'linewidth', 2);
plot(max(Wdata) * [0.1 0.1], max(EE) * [0 1], 'k--');
plot(max(Wdata) * [midpoint midpoint], max(EE) * [0 1], 'k--');
plot(max(Wdata) * [0.9 0.9], max(EE) * [0 1], 'k--');
text(max(Wdata)*0.2, 200, 'Low W')
text(max(Wdata)*0.6, 200, 'High W')

text(max(Wdata)*0.2, 1400, sprintf('y=%2.1fx+%2.1f', mdl1k, mdl1m))
text(max(Wdata)*0.6, 1400, sprintf('y=%2.1fx+%2.1f', mdl2k, mdl2m))

ylim([0 inf])
xlim([0 inf])


mdl2k/mdl1k

