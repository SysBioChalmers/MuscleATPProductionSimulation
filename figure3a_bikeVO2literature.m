close all
%Based on published data from:
%https://www.ncbi.nlm.nih.gov/pubmed/12510865

subjectData = [
    1	24	79	184	23.33	42.3
    2	22	80	182	24.15	53.7
    3	22	74	180	22.84	53.2
    4	22	73	176	23.57	50.1
    5	22	72	183	21.50	50.0
    6	21	59	171	20.18	42.8
    7	24	71	180	21.94	46.8
    8	25	73	172	24.68	51.8
    9	22	70	180	21.60	56.3
    10	26	60	174	19.82	58.8
    11	24	77	179	24.03	44.7
    12	26	86	186	24.86	49.4
    13	28	71	176	22.92	55.4
    14	21	65	173	21.72	55.8
    15	24	68	177	21.71	53.6
    16	22	70	183	20.90	50.0
    17	23	79	176	25.50	53.1
    18	21	73	188	20.65	52.2
    19	29	68	176	21.95	59.4
    20	22	77	180	23.77	48.4
    21	29	87	180	26.85	40.6
    ];

%
muscleType =[
1	53.9	46.1	33.1	13.0
2	64.9	35.1	35.1	0
3	55.3	44.7	44.7	0
4	52.7	47.3	41.8	5.5
5	38.7	61.3	47.2	14.1
6	61.1	38.9	38.9	0
7	49.2	50.8	41.0	9.8
8	82.5	17.5	17.5	0
9	68.3	31.7	31.7	0
10	66.0	34.0	34.0	0
11	57.7	42.3	38.3	4.0
12	16.1	83.9	54.3	29.6
13	44.2	55.8	34.7	21.1
14	45.6	54.4	45.9	8.5
15	46.0	54.0	46.9	7.1
16	57.9	42.1	29.4	12.7
17	70.2	29.8	29.8	0
18	43.3	56.7	45.6	11.1
19	79.2	20.8	20.8	0
20	71.7	28.3	28.3	0
21	69.6	30.4	30.4	0
];

%subject W %ofvO2max W max
LT = [1	150	59.5	250
2	270	75.1	345
3	210	67.5	310
4	210	70.4	337
5	120	46.5	265
6	120	56.4	240
7	120	49.4	258
8	180	60.5	305
9	150	51.9	303
10	150	54.0	280
11	150	59.7	273
12	150	48.4	317
13	150	53.8	292
14	150	54.2	282
15	150	53.0	270
16	150	59.2	270
17	150	47.7	330
18	120	42.8	302
19	210	63.0	330
20	180	63.4	305
21	180	67.8	292
];

% subject intercept slope
linFit = [1	451.6	10.227
2	552.86	9.81
3	426.71	10.394
4	710.29	8.62
5	715.5	8.137
6	335	8.803
7	462.5	9.83
8	392.6	10.556
9	521.4	9.7
10	454.2	9.413
11	552.8	9.94
12	612.6	9.233
13	404.8	11.12
14	454.5	10.103
15	577.8	8.827
16	485.2	10.487
17	633.3	8.743
18	629.5	8.497
19	179.14	11.186
20	577.73	9.668
21	405.33	10.562
];

observedVsExpected = [
    1	3344	3008	336
    2	4299	3987	312
    3	3933	3649	284
    4	3660	3615	45
    5	3600	2872	728
    6	2526	2448	78
    7	3324	2999	325
    8	3779	3612	167
    9	3944	3461	483
    10	3525	3090	435
    11	3438	3266	172
    12	4252	3539	713
    13	3931	3652	279
    14	3626	3304	322
    15	3643	2961	682
    16	3503	3317	186
    17	4197	3519	678
    18	3811	3196	615
    19	4038	3871	167
    20	3726	3526	200
    21	3530	3489	41
    ];


slope1 = linFit(:,3);
deltaW = LT(:,4)-LT(:,2);
O2LT = linFit(:,2) + LT(:,2) .* linFit(:,3);
O2VO2max = observedVsExpected(:,2);
slope2 = (O2VO2max-O2LT)./deltaW;

deltas = [slope1;slope2];
groups = [ones(length(slope1),1); 2 * ones(length(slope2),1)];

groupNames = {'low W', 'high W'};
h = boxplot(deltas, groups, 'Orientation',  'vertical',  'Colors', 'k', 'Widths', 1);
set(h(7,:),'Visible','off')
set(gca,'xticklabel', groupNames)  
ylim([0 17])
ylabel('dO2/dW')

color = [67 116 160; 151 185 224]/255;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    currentColor = color((mod(j,2)+1),:);
    a = patch(get(h(j),'XData'),get(h(j),'YData'), currentColor);
    uistack(a,'bottom');
end
%uistack(a,'bottom');
p = ranksum(deltas(groups==1),deltas(groups==2));
%[h p] = ttest(deltas(groups==1),deltas(groups==2));

text(1, 4, sprintf('p=%2.2e', p))
ylim([0 16])
color2 = [215 86 40]/256;
hold all
rng('default') % set random generator to same value each time
xvals = 0.5*rand(length(slope1),1);
xvals = xvals-mean(xvals);

% for i = 1:length(xvals)
%     tmp = plot([1 2] + xvals(i), [slope1(i) slope2(i)], '-', 'color', color2);
%     tmp.Color(4) = 0.5;
% end

scatter(1 + xvals, slope1, 20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
scatter(2 + xvals, slope2, 20, 'MarkerFaceColor', color2,'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)


text(0.3, median(slope1), sprintf('%2.1f', median(slope1)))
text(1.3, median(slope2), sprintf('%2.1f', median(slope2)))

%%
addpath('src2')
Wvalues = [median(slope1) median(slope2)];
Ovalue = 2 * mlToMol(Wvalues)/3600;
ATPvalue = (1/27)/1000;
ATPvalue./Ovalue




