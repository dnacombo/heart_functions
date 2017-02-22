%% Setup

% add the folder with heart_SV_calc.m and heart_peak_detect.m to the path
% (you need to edit this)
add_fieldtrip
addpath(cdhome('heart'))

file = 'MERGE_CUT_pilote1_damier_EBI_12_5_kHz_tsss.fif';

% thoracic length (cm)
L = 30;

%% Read the data using fieldtrip
% Simply read a dataset with a cardiac data.

cfg                         = [];
cfg.dataset                 = file;
cfg.channel                 = {'EBI_-_Magnitude' 'ECG1' 'HeartSound'};
data_raw                    = ft_preprocessing(cfg);

%% Detecting heart peaks
% If there is an ECG channel, we can detect heart peaks P Q R S T.
% Everything is based on detecting the R peaks.

cfg = [];
cfg.channel = 'ECG1';

[HeartBeats] = heart_peak_detect(cfg,data_raw);

HeartBeats
% retrieve R peaks with [HeartBeats.R_sample] or [HeartBeats.R_time]

%% Using visual feedback to improve detection
%

cfg = [];
cfg.channel = 'ECG1';
cfg.plotall = 'yes';
[HeartBeats] = heart_peak_detect(cfg,data_raw);


%% Stroke volume
% Simplest call without any feedback
cfg = [];
cfg.interactive = 'no';
cfg.L = L;

[HeartBeats] = heart_SV_calc(cfg,data_raw);
% retrieve Stroke volume with [HeartBeats.SV]

%% A simple plot
% A histogram of stroke volumes extracted at each heart beat

% Stroke volume in mL
StrokeVol = [HeartBeats.SV];

figure;
hist(StrokeVol,20)
xlabel('Stroke volume (mL)')
ylabel('Number of strokes');

%% Use feedback to remove outliers

cfg = [];
cfg.L = L;

[HeartBeats] = heart_SV_calc(cfg,data_raw);

%% A new simple plot
% A histogram of stroke volumes extracted at each heart beat

StrokeVol = [HeartBeats.SV];

figure;
hist(StrokeVol,20)
xlabel('Stroke volume (mL)')
ylabel('Number of strokes');
%% Compute other useful values
% Heart rate, IBI, Cardiac output...

% R peaks in seconds 
R = [HeartBeats.R_time];

% Average heart rate in pulse/min
HR = numel(R) / (R(end) - R(1)) * 60;

% IBI in seconds
IBI = diff(R);

% Instantaneous Heart Rate at each HB
IHR = 1./IBI.*60;

% Cardiac output in L/min
CO = mean(StrokeVol) * HR / 1000;

%% Without preloading the data

cfg = [];
cfg.dataset = file;
cfg.L = L;

[HeartBeats] = heart_SV_calc(cfg);



