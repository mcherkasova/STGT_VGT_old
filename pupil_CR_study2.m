clc;clear; close all;

% in directories
% change 'im' to 'gt' or 'st' in file names and directories, depending on
% which group is being analyzed
pupil_dirName = '***\Pupil\study2\pupil\im';
marker_dirName = '***\Pupil\study2\markers\im';
% out directories
timecourses_dirName = '***\Pupil\study2\timecourses';

files = dir(fullfile(pupil_dirName, '*.mat'));
files = {files.name};

% AUC file
fID = fopen('***\Pupil\study2\CS_AUC_im.txt','a');
fprintf(fID,'%12s %12s %12s %12s %12s %12s\r\n',...
    'ID','CS','Hemiblock','trial','AUC','AUC_500lag');

for z=1:numel(files)
%% load 
pupil_filename = fullfile(pupil_dirName,files{z});
marker_filename = fullfile(marker_dirName, files{z});
load(pupil_filename);
load(marker_filename);
id=str2num(pupil_filename(end-6:end-4));
%% Filter pupil
pupil_unfiltered = pupil;
hf = 4*2/1000;
[b,a] = butter(2, hf);
pupil = filter(b,a,pupil_unfiltered);

%% Define onsets
%Hemiblock 1
%CS plus
bl_plus_h1_start = onsets{1};
bl_plus_h1_end = onsets{5};
CSplus_h1_start = onsets{9};
CSplus_h1_lag = onsets{13};
FBplus_h1_start = onsets{17};
FBplus_h1_end = FBplus_h1_start+2000;

%CS minus
bl_minus_h1_start = onsets{3};
bl_minus_h1_end = onsets{7};
CSminus_h1_start = onsets{11};
CSminus_h1_lag = onsets{15};
FBminus_h1_start = onsets{19};
FBminus_h1_end = FBminus_h1_start+2000;

%Hemiblock 2
%CS plus
bl_plus_h2_start = onsets{2};
bl_plus_h2_end = onsets{6};
CSplus_h2_start = onsets{10};
CSplus_h2_lag = onsets{14};
FBplus_h2_start = onsets{18};
FBplus_h2_end = FBplus_h2_start+2000;

%CS minus
bl_minus_h2_start = onsets{4};
bl_minus_h2_end = onsets{8};
CSminus_h2_start = onsets{12};
CSminus_h2_lag = onsets{16};
FBminus_h2_start = onsets{20};
FBminus_h2_end = FBminus_h2_start+2000;

%% Define per trial timecourses and calulate AUCs and mean timecourses
%Full trial
%Find minimum duration
for kk=1:10
    dur1=FBplus_h1_end(kk)-bl_plus_h1_start(kk);
    dur2=FBminus_h1_end(kk)-bl_minus_h1_start(kk);
    dur3=FBplus_h2_end(kk)-bl_plus_h2_start(kk);
    dur4=FBminus_h2_end(kk)-bl_minus_h2_start(kk);
    durarr=[dur1, dur2, dur3,dur4];
    mindur = min(durarr);
    durs(kk)= mindur;
end

mindur = min(durs);

[FullCSplusH1_timecourse_mat, FullCSplusH1_mat_pos, FullCSplusH1_AUCs,FullCSplusH1_timecourse_mean] = parsepupil(pupil,...
    bl_plus_h1_start,FBplus_h1_end,mindur);
[FullCSminusH1_timecourse_mat, FullCSminusH1_mat_pos, FullCSminusH1_AUCs,FullCSminusH1_timecourse_mean] = parsepupil(pupil,...
    bl_minus_h1_start,FBminus_h1_end,mindur);
[FullCSplusH2_timecourse_mat, FullCSplusH2_mat_pos, FullCSplusH2_AUCs,FullCSplusH2_timecourse_mean] = parsepupil(pupil,...
    bl_plus_h2_start,FBplus_h2_end,mindur);
[FullCSminusH2_timecourse_mat, FullCSminusH2_mat_pos, FullCSminusH2_AUCs,FullCSminusH2_timecourse_mean] = parsepupil(pupil,...
    bl_minus_h2_start,FBminus_h2_end,mindur);

FullCSplusH1timecourse(z,:)=[id,{FullCSplusH1_timecourse_mean}]; % Compile mean timecourse matrices for full trial
FullCSminusH1timecourse(z,:)=[id,{FullCSminusH1_timecourse_mean}];
FullCSplusH2timecourse(z,:)=[id,{FullCSplusH2_timecourse_mean}];
FullCSminusH2timecourse(z,:)=[id,{FullCSminusH2_timecourse_mean}];

%CS
[CSplusH1_timecourse_mat, CSplusH1_mat_pos, CSplusH1_AUCs,CSplusH1_timecourse_mean] = parsepupil(pupil,...
    CSplus_h1_start,FBplus_h1_start,6900);
[CSminusH1_timecourse_mat, CSminusH1_mat_pos, CSminusH1_AUCs,CSminusH1_timecourse_mean] = parsepupil(pupil,...
    CSminus_h1_start,FBminus_h1_start,6900);
[CSplusH2_timecourse_mat, CSplusH2_mat_pos, CSplusH2_AUCs,CSplusH2_timecourse_mean] = parsepupil(pupil,...
    CSplus_h2_start,FBplus_h2_start,6900);
[CSminusH2_timecourse_mat, CSminusH2_mat_pos, CSminusH2_AUCs,CSminusH2_timecourse_mean] = parsepupil(pupil,...
    CSminus_h2_start,FBminus_h2_start,6900);

%CS with a 500ms lag
[CSplusH1lag_timecourse_mat, CSplusH1lag_mat_pos, CSplusH1_lag_AUCs,CSplusH1lag_timecourse_mean] = parsepupil(pupil,...
    CSplus_h1_lag,FBplus_h1_start,6400);
[CSminusH1lag_timecourse_mat, CSminusH1lag_mat_pos, CSminusH1_lag_AUCs,CSminusH1lag_timecourse_mean] = parsepupil(pupil,...
    CSminus_h1_lag,FBminus_h1_start,6400);
[CSplusH2lag_timecourse_mat, CSplusH2lag_mat_pos, CSplusH2_lag_AUCs,CSplusH2lag_timecourse_mean] = parsepupil(pupil,...
    CSplus_h2_lag,FBplus_h2_start,6400);
[CSminusH2lag_timecourse_mat, CSminusH2lag_mat_pos, CSminusH2_lag_AUCs,CSminusH2lag_timecourse_mean] = parsepupil(pupil,...
    CSminus_h2_lag,FBminus_h2_start,6400);

id_array = repmat(id,40,1); % Create matrix for AUCs to be written into text file
cs_array=[repmat(1,10,1);repmat(-1,10,1); repmat(1,10,1);repmat(-1,10,1)];
hemiblock_array=[repmat(1,20,1); repmat(2,20,1)];
trail_array = [onsets{21};onsets{23};onsets{22};onsets{24}];
AUC_mat = [CSplusH1_AUCs'; CSminusH1_AUCs'; CSplusH2_AUCs'; CSminusH2_AUCs']; 
AUC_mat_lag = [CSplusH1_lag_AUCs';CSminusH1_lag_AUCs';CSplusH2_lag_AUCs';CSminusH2_lag_AUCs']; 
AUC_mat = [id_array, cs_array, hemiblock_array, trail_array, AUC_mat, AUC_mat_lag];
fprintf(fID,'%4d %4d %4d %4d %4f %4f\r\n',AUC_mat'); %Write to txt file

CSplusH1timecourse(z,:)=[id,CSplusH1_timecourse_mean]; % Compile mean timecourse matrices for CS epoch
CSminusH1timecourse(z,:)=[id,CSminusH1_timecourse_mean];
CSplusH2timecourse(z,:)=[id,CSplusH2_timecourse_mean];
CSminusH2timecourse(z,:)=[id,CSminusH2_timecourse_mean];

% Feedback
[FBplusH1_timecourse_mat, FBplusH1_mat_pos, FBplusH1_AUCs,FBplusH1_timecourse_mean] = parsepupil(pupil,...
    FBplus_h1_start,FBplus_h1_end,1000);
[FBminusH1_timecourse_mat, FBminusH1_mat_pos, FBminusH1_AUCs,FBminusH1_timecourse_mean] = parsepupil(pupil,...
    FBminus_h1_start,FBminus_h1_end,1000);
[FBplusH2_timecourse_mat, FBplusH2_mat_pos, FBplusH2_AUCs,FBplusH2_timecourse_mean] = parsepupil(pupil,...
    FBplus_h2_start,FBplus_h2_end,1000);
[FBminusH2_timecourse_mat, FBminusH2_mat_pos, FBminusH2_AUCs,FBminusH2_timecourse_mean] = parsepupil(pupil,...
    FBminus_h2_start,FBminus_h2_end,1000);

FBplusH1timecourse(z,:)=[id,FBplusH1_timecourse_mean]; % Compile mean timecourse matrices for Feedback epoch
FBminusH1timecourse(z,:)=[id,FBminusH1_timecourse_mean];
FBplusH2timecourse(z,:)=[id,FBplusH2_timecourse_mean];
FBminusH2timecourse(z,:)=[id,FBminusH2_timecourse_mean];

clearvars -except fID files pupil_dirName marker_dirName timecourses_dirName...
     FullCSplusH1timecourse FullCSminusH1timecourse FullCSplusH2timecourse FullCSminusH2timecourse...
     CSplusH1timecourse CSminusH1timecourse CSplusH2timecourse CSminusH2timecourse...
     FBplusH1timecourse FBminusH1timecourse FBplusH2timecourse FBminusH2timecourse
    
end

save('***\Pupil\study2\timecourses\CSminusH1timecourse_im.mat','CSminusH1timecourse');
save('***\Pupil\study2\timecourses\CSminusH2timecourse_im.mat','CSminusH2timecourse');
save('***\Pupil\study2\timecourses\CSplusH1timecourse_im.mat','CSplusH1timecourse');
save('***\Pupil\study2\timecourses\CSplusH2timecourse_im.mat','CSplusH2timecourse');
save('***\Pupil\study2\timecourses\FBminusH1timecourse_im.mat','FBminusH1timecourse');
save('***\Pupil\study2\timecourses\FBminusH2timecourse_im.mat','FBminusH2timecourse');
save('***\Pupil\study2\timecourses\FBplusH1timecourse_im.mat','FBplusH1timecourse');
save('***\Pupil\study2\timecourses\FBplusH2timecourse_im.mat','FBplusH2timecourse');
save('***\Pupil\study2\timecourses\FullCSminusH1timecourse_im.mat','FullCSminusH1timecourse');
save('***\Pupil\study2\timecourses\FullCSminusH2timecourse_im.mat','FullCSminusH2timecourse');
save('***\Pupil\study2\timecourses\FullCSplusH1timecourse_im.mat','FullCSplusH1timecourse');
save('***\Pupil\study2\timecourses\FullCSplusH2timecourse_im.mat','FullCSplusH2timecourse');

fclose('all');


function [timecourse_mat, timecourse_mat_pos, AUCs,timecourse_mean] = parsepupil(pupil,start_marker,...
    end_marker,timecourse_length)
% Creates matrices with trial-wise pupil timecourses delineated by start
% and end markers. Also produces an area under the curve (AUC) array for
% these timecourses using positive values (by adding a constant 10 to
% timecourses). In addition, produces a mean timecourse for the 10 trials.

for j=1:10
    tr = start_marker(j):end_marker(j);
    pupil_per_trial = pupil(tr);
    pupil_per_trial = pupil_per_trial(1:timecourse_length);
    pupil_per_trial_pos = pupil_per_trial+10;
    AUC = trapz(pupil_per_trial_pos);
    timecourse_mat(j,:) = pupil_per_trial;
    timecourse_mat_pos(j,:) = pupil_per_trial_pos;
    AUCs(j)= AUC;
end

timecourse_mean = nanmean(timecourse_mat,1);% mean pupil size
end


