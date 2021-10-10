clc;clear; close all;

%% Calculate means and SEs
%CS phase
mean_CSminusH1 = nanmean(CSminusH1timecourse(:,2:end),1);
se_CSminusH1 = nanstd(CSminusH1timecourse(:,2:end),1)/sqrt(19);

mean_CSplusH1 = nanmean(CSplusH1timecourse(:,2:end),1);
se_CSplusH1 = nanstd(CSplusH1timecourse(:,2:end),1)/sqrt(19);

mean_CSminusH2 = nanmean(CSminusH2timecourse(:,2:end),1);
se_CSminusH2 = nanstd(CSminusH2timecourse(:,2:end),1)/sqrt(19);

mean_CSplusH2 = nanmean(CSplusH2timecourse(:,2:end),1);
se_CSplusH2 = nanstd(CSplusH2timecourse(:,2:end),1)/sqrt(19);

%FB phase
mean_FBminusH1 = nanmean(FBminusH1timecourse(:,2:end),1);
se_FBminusH1 = nanstd(FBminusH1timecourse(:,2:end),1)/sqrt(19);

mean_FBplusH1 = nanmean(FBplusH1timecourse(:,2:end),1);
se_FBplusH1 = nanstd(FBplusH1timecourse(:,2:end),1)/sqrt(19);

mean_FBminusH2 = nanmean(FBminusH2timecourse(:,2:end),1);
se_FBminusH2 = nanstd(FBminusH2timecourse(:,2:end),1)/sqrt(19);

mean_FBplusH2 = nanmean(FBplusH2timecourse(:,2:end),1);
se_FBplusH2 = nanstd(FBplusH2timecourse(:,2:end),1)/sqrt(19);

% CS plus feedback
response=nan(1,5);
meanCSminus_FB_H1 = [mean_CSminusH1,response,mean_FBminusH1];
seCSminus_FB_H1 = [se_CSminusH1,response,se_FBminusH1];
meanCSplus_FB_H1 = [mean_CSplusH1,response,mean_FBplusH1];
seCSplus_FB_H1 = [se_CSplusH1,response,se_FBplusH1];

meanCSminus_FB_H2 = [mean_CSminusH2,response,mean_FBminusH2];
seCSminus_FB_H2 = [se_CSminusH2,response,se_FBminusH2];
meanCSplus_FB_H2 = [mean_CSplusH2,response,mean_FBplusH2];
seCSplus_FB_H2 = [se_CSplusH2,response,se_FBplusH2];

%% Plots
% CS phase
x=1:6901;%adjust depending on the study
figure(1)
[hl1, hp1] = boundedline(x, mean_CSminusH1, se_CSminusH1, 'k', x, mean_CSplusH1, se_CSplusH1, 'b', 'transparency', .1, 'alpha'); %plot
%outlinebounds(hl1,hp1);
title('GT Hemiblock 1')
xlabel('Time (ms)')
ylabel('Pupil Size (z scores)')
legend('CS-', 'CS+')

figure(2)
[hl1, hp1] = boundedline(x, mean_CSminusH2, se_CSminusH2, 'k', x, mean_CSplusH2, se_CSplusH2, 'b', 'transparency', .1, 'alpha'); %plot
%outlinebounds(hl1,hp1);
title('GT Hemiblock 2')
xlabel('Time (ms)')
ylabel('Pupil Size (z scores)')
legend('CS-', 'CS+')

% Feedback phase
x=1:1000;%adjust depending on the study
figure(1)
[hl1, hp1] = boundedline(x, mean_FBminusH1, se_FBminusH1, 'k', x, mean_FBplusH1, se_FBplusH1, 'b', 'transparency', .1, 'alpha'); %plot
%outlinebounds(hl1,hp1);
title('Feedback GT Hemiblock 1')
xlabel('Time (ms)')
ylabel('Pupil Size (z scores)')
legend('CS-', 'CS+')

figure(2)
[hl1, hp1] = boundedline(x, mean_FBminusH2, se_FBminusH2, 'k', x, mean_FBplusH2, se_FBplusH2, 'b', 'transparency', .1, 'alpha'); %plot
%outlinebounds(hl1,hp1);
title('Feedback GT Hemiblock 2')
xlabel('Time (ms)')
ylabel('Pupil Size (z scores)')
legend('CS-', 'CS+')

% Combined CS and feedback phase
x=1:7905;%adjust depending on the study
ax(1) = subplot(2,2,1);
[l, p] = boundedline(x, meanCSminus_FB_H1, seCSminus_FB_H1, 'k', x, meanCSplus_FB_H1, seCSplus_FB_H1, 'b', 'transparency', .1, 'alpha'); %plot
%outlinebounds(hl1,hp1);
title('ST Hemiblock 1')
xlabel('Time (ms)')
ylabel('Pupil Size (z scores)')
legend('CS-', 'CS+', 'Location','southeast')
ylim([-1.5 1])

ax(1) = subplot(2,2,2);
boundedline(x, meanCSminus_FB_H2, seCSminus_FB_H2, 'k', x, meanCSplus_FB_H2, seCSplus_FB_H2, 'b', 'transparency', .1, 'alpha'); %plot
%outlinebounds(hl1,hp1);
title('ST Hemiblock 2')
xlabel('Time (ms)')
ylabel('Pupil Size (z scores)')
%legend('CS-', 'CS+', 'Location','south')
ylim([-1.5 1])