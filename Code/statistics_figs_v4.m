% Paper Figures 1-5
function [] = statistics_figs_v4(LAB_SW,IMU_all,VICON,ClosestHome_WS,HOME_WS,HOME_WS_date,HOME_assist,lme_ast)

%% set-up for lab data
controls_SWdata_all = LAB_SW(10:15);
controls_imudata_all = IMU_all(10:15);
controls_vicondata_all = VICON;
patient_SWdata_all = LAB_SW(1:9);
patient_imudata_all = IMU_all(1:9);

% combine controls and patients lab imu data
labimudata_all = [patient_imudata_all,controls_imudata_all];

%% set-up for home data
homedata_all = ClosestHome_WS.val;

%% Lab Comparisons Boxplot - Controls (Figure 2)
labdata_controls = [];
labdata_groupC = [];
j = 1; % counter for IMU data
k = 2; % counter for SW data
m = 3; % counter for Vicon data
n = 4; % counter for nan data
% add a grouping of nans in between sets
space = nan(5,1);
for i = 1:length(controls_imudata_all)
    labdata_controls = [labdata_controls;controls_imudata_all{i}';controls_SWdata_all{i}';controls_vicondata_all{i}';space];
    labdata_groupC = [labdata_groupC;repmat(j,length(controls_imudata_all{i}),1);repmat(k,length(controls_SWdata_all{i}),1);repmat(m,length(controls_vicondata_all{i}),1);repmat(n,length(space),1)];
    j = j + 4;
    k = k + 4;
    m = m + 4;
    n = n + 4;
end
colors = [0 0 0; .7 0 0; 0 0 .7; 0 0 0];
color_all = repmat(colors,length(controls_imudata_all),1);
lbls_con = repmat({'M','S','V',''},1,length(controls_imudata_all));
h = figure('Position',[0 0 300 300]);
boxplot(labdata_controls,labdata_groupC,'labels',lbls_con,'LabelOrientation','inline','Colors',color_all,'widths',.8)
o = findobj(gcf,'tag','Outlier');
set(o,'Color','k')
ylabel('Walking Speed (m/s)')
ax = gca;
ax.FontSize = 12;

% add second level labels
lowerLevel = 0.925;
lowerLabels = {'C1','C2','C3','C4','C5','C6'};
text([2, 6, 10, 14, 18, 22], repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','FontSize',12)
if ~(exist('paperfigures')==7)
    mkdir('paperfigures')
end
savepdffigure(h,'paperfigures/FIG2_ControlHomeLabComparison')

%% Lab Comparisons Boxplot -Patients (Figure 1)
labdata_patients = [];
labdata_group = [];
j = 1; % counter for IMU data
k = 2; % counter for SW data
n = 3; % counter for nan data
% add a grouping of nans in between sets
space = nan(5,1);
for i = 1:length(patient_imudata_all)
    labdata_patients = [labdata_patients;patient_imudata_all{i}';patient_SWdata_all{i}';space];
    labdata_group = [labdata_group;repmat(j,length(patient_imudata_all{i}),1);repmat(k,length(patient_SWdata_all{i}),1);repmat(n,length(space),1)];
    j = j + 3;
    k = k + 3;
    n = n + 3;
end

colors = [0 0 0; .7 0 0; 0 0 0];
color_all = repmat(colors,length(patient_imudata_all),1);
lbls_pat = repmat({'M','S',''},1,length(patient_imudata_all));
h = figure('Position',[0 0 300 300]);
boxplot(labdata_patients,labdata_group,'labels',lbls_pat,'LabelOrientation','horizontal','Colors',color_all,'widths',.8)
ylabel('Walking Speed (m/s)')
ax = gca;
ax.FontSize = 12;

% add second level labels
lowerLevel = 0.25;
lowerLabels = {'P3','P5','P7','P8','P9','P10','P11','P12','P13'};
text([1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5, 25.5], repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
    'VerticalAlignment','Top','HorizontalAlignment','Center','FontSize',12)
savepdffigure(h,'paperfigures/FIG1_PatientHomeLabComparison')

%% Scatter plot of mean clinic and home walking speed for all subjects (Figure 3)
% find mean values
clinM = zeros(length(labimudata_all),1);
homeM = zeros(length(homedata_all),1);
for i = 1:length(labimudata_all)
    clinM(i) = mean(labimudata_all{i});
    homeM(i) = mean(homedata_all{i});
end

h = figure('Position',[0 0 250 250]);
% patients
scatter(homeM(1:9),clinM(1:9),[],'ko');
hold on
scatter(homeM(10:15),clinM(10:15),[],'k+');
lndiag = linspace(0,1.5,100);
plot(lndiag,lndiag,':k');
xlim([0 1.5])
ylim([0 1.5])

ylabel('Lab IMU Walking Speed (m/s)','FontSize',12);
xlabel('Home IMU Walking Speed (m/s)','FontSize',12);
legend('Patient','Control','FontSize',12,'Location','SouthEast');
axis square

savepdffigure(h,'paperfigures/FIG3_LabHomeScatter');

%% Walking speed over time for patients (Figure 4)
h=figure('Position',[100 100 1000 500]);
clear ws_mean days_mean assist_change
for i = 1:9
    sp = subplot(3,3,i);
    ws_mean{i} = cellfun(@mean,HOME_WS(i,:));
    ws_mean{i}(isnan(ws_mean{i})) = [];
    days_mean{i} = [HOME_WS_date{i,:}];
    data_tbl = table(ws_mean{i}',days_mean{i}'-min(days_mean{i}),'VariableNames',{'WalkingSpeed','Days'});
    lme_ws{i} = fitlm(data_tbl,'WalkingSpeed ~ Days');

    if strcmp(ClosestHome_WS.code{i},HOME_assist{1,1}) %SC07
        [a,b,c] = unique(HOME_assist{1,2},'rows','stable');
        assist_change{1} = HOME_assist{1,4}(b);
    end
    if strcmp(ClosestHome_WS.code{i},unique(HOME_assist{2,1})) %SC08
        [a,b,c] = unique(HOME_assist{2,2},'rows','stable');
        assist_change{2} = HOME_assist{2,4}(b);
    end
    if strcmp(ClosestHome_WS.code{i},unique(HOME_assist{3,1})) %SC11
        [a,b,c] = unique(HOME_assist{3,2},'rows','stable');
        assist_change{3} = HOME_assist{3,4}(b);
    end

    scatter(days_mean{i}-min(days_mean{i}),ws_mean{i},10);
    hold on
    ln_mean = lme_ws{i}.Coefficients.Estimate(2)*(days_mean{i}-min(days_mean{i})) +...
        lme_ws{i}.Coefficients.Estimate(1);
    if   lme_ws{i}.Coefficients.pValue < 0.05 % significant change over time
        plot(days_mean{i}-min(days_mean{i}),ln_mean,'Color','r','LineWidth',1,'LineStyle','-')
    end
    sp.FontSize = 12;
    ylim([0.3 1.6])
    xlim([0 200])
    xlabel('Days','FontSize',12)
    ylabel('Walking Speed (m/s)','FontSize',12)

    title([lowerLabels{i} ': $v=' num2str(round(lme_ws{i}.Coefficients.Estimate(2),3,'significant'))...
        'x_D+' num2str(round(  lme_ws{i}.Coefficients.Estimate(1),3,'significant')) '$'],...
        'FontSize',12,'Interpreter','latex')
end
savepdffigure(h,'paperfigures/FIG4_WSOverTime');

%% Mobility assistance change (Figure 5)

% colors
color = parula;
color2 = winter;
blue = [color2(1,:);color(70,:)];
color3 = hot;
orange = [color3(145,:);0.8500, 0.3250, 0.0980];

%SC07
ws7m = ws_mean{3};
days7m = days_mean{3}-min(days_mean{3});
ad7m = assist_change{1};
ws7m_noad = ws7m(ad7m==0);
ws7m_ad = ws7m(ad7m>0);
days7m_noad = days7m(ad7m==0);
days7m_ad = days7m(ad7m>0);
h=figure
subplot(3,1,1)
scatter(days7m_noad,ws7m_noad,[],blue(2,:),'o')
hold on
ln7m_noad = lme_ast{1}.Coefficients.Estimate(1) + ...
    lme_ast{1}.Coefficients.Estimate(2)*days7m_noad;
plot(days7m_noad,ln7m_noad,'--','LineWidth',1,'Color',blue(1,:))
scatter(days7m_ad,ws7m_ad,[],orange(1,:),'x')
ln7m_ad = (lme_ast{1}.Coefficients.Estimate(1) + lme_ast{1}.Coefficients.Estimate(3)) + ...
    (lme_ast{1}.Coefficients.Estimate(2) + lme_ast{1}.Coefficients.Estimate(4))*...
    days7m_ad;
plot(days7m_ad,ln7m_ad,'--','LineWidth',1,'Color',orange(2,:))
ylim([0.3 1.6])
ylabel('Walking speed (m/s)')
xlim([0 190])
title('P7')
legend('Walking speed without assist device','Change in walking speed without assist device','Walking speed with assist device','Change in walking speed with assist device')

%SC08
subplot(3,1,2)
ws8m = ws_mean{4};
days8m = days_mean{4}-min(days_mean{4});
ad8m = assist_change{2};
ws8m_noad = ws8m(ad8m==0);
ws8m_ad = ws8m(ad8m>0);
days8m_noad = days8m(ad8m==0);
days8m_ad = days8m(ad8m>0);
scatter(days8m_noad,ws8m_noad,[],blue(2,:),'o')
hold on
ln8m_noad = lme_ast{2}.Coefficients.Estimate(1) + ...
    lme_ast{2}.Coefficients.Estimate(2)*days8m_noad;
plot(days8m_noad,ln8m_noad,'--','LineWidth',1,'Color',blue(1,:))
scatter(days8m_ad,ws8m_ad,[],orange(1,:),'x')
ln8m_ad = (lme_ast{2}.Coefficients.Estimate(1) + lme_ast{2}.Coefficients.Estimate(3)) + ...
    (lme_ast{2}.Coefficients.Estimate(2) + lme_ast{2}.Coefficients.Estimate(4))*...
    days8m_ad;
plot(days8m_ad,ln8m_ad,'--','LineWidth',1,'Color',orange(2,:))
ylim([0.3 1.6])
ylabel('Walking speed (m/s)')
xlim([0 190])
title('P8')

%SC11
subplot(3,1,3)
ws11m = ws_mean{7};
days11m = days_mean{7}-min(days_mean{7})
ad11m = assist_change{3};
ws11m_noad = ws11m(ad11m<2);
ws11m_ad = ws11m(ad11m==2);
days11m_noad = days11m(ad11m<2);
days11m_ad = days11m(ad11m==2);
scatter(days11m_noad,ws11m_noad,[],blue(2,:),'o')
hold on
ln11m_noad = lme_ast{3}.Coefficients.Estimate(1) + ...
    lme_ast{3}.Coefficients.Estimate(2)*days11m_noad;
plot(days11m_noad,ln11m_noad,'--','LineWidth',1,'Color',blue(1,:))
scatter(days11m_ad,ws11m_ad,[],orange(1,:),'x')
ln11m_ad = (lme_ast{3}.Coefficients.Estimate(1) + lme_ast{3}.Coefficients.Estimate(3)) + ...
    (lme_ast{3}.Coefficients.Estimate(2) + lme_ast{3}.Coefficients.Estimate(4))*...
    days11m_ad;
plot(days11m_ad,ln11m_ad,'--','LineWidth',1,'Color',orange(2,:))
ylim([0.3 1.6])
xlabel('Days after 1st home recording')
ylabel('Walking speed (m/s)')
xlim([0 190])
title('P11')
savepdffigure(h,'paperfigures/FIG5_WalkAssist');


