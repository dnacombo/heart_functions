function HeartBeats = heart_SV_calc(cfg,data)

% HeartBeats = heart_SV_calc(cfg,data)
%
%

def             = [];
def.ECG.channel = 'ECG1';
def.ECG.magfactor = -4000; % just for display along with Z0
def.ECG.show    = 'yes';

def.Rpeakcfg    = [];

def.Z0.channel  = 'EBI_-_Magnitude';
def.Z0.lpfilter = 'yes';
def.Z0.lpfreq   = 10;
def.Z0.lpfilttype = 'but';
def.Z0.factor   = 5;
def.Z0.show     = 'yes';

def.sound.channel = 'HeartSound';
def.sound.bpfilter = 'yes';
def.sound.bpfreq = [40 60];
def.sound.tubelength = 152;

def.LVE_startfracmax_dZdt = 0.15;
def.rho         = 135;
def.L           = 30;
def.Z           = 40;
def.SV          = @(rho,L,Z,LVET,dZdt_max_after_HB) rho .* L .^ 2 ./ Z .^ 2 .* LVET .* dZdt_max_after_HB;

def.plotallmeas = 'no';
def.interactive = 'yes';

cfg = setdef(cfg,def);

hasdata = exist('data','var');

if ~hasdata % try to read data from disk using cfg.dataset
    data = ft_preprocessing(cfg);
end
time = data.time{1};

iECG = ft_channelselection(cfg.ECG.channel,data.label);
iZ0 = ft_channelselection(cfg.Z0.channel,data.label);
isound = ft_channelselection(cfg.sound.channel,data.label);

if isempty(iECG)
    error('no ECG channel found')
end
if numel(iECG) > 1
    disp(['ECG channel is average of ' sprintf('%s ',iECG{:})])
end
ECG = mean(data.trial{1}(chnb(iECG,data.label),:),1);

if isempty(iZ0)
    error('no impedance channel found')
end
if numel(iZ0) > 1
    error('I don''t know how to work with several impedance channels')
end

if isempty(isound)
    disp('No heart sound channel provided')
    hassound = 0;
else
    hassound = 1;
end
if numel(isound) > 1
    error('I don''t know how to work with several heart sound channels')
end

cfg.Rpeakcfg.fsample = data.fsample;
[HB, HB_time] = heart_peak_detection(ECG,cfg.Rpeakcfg);
% remove last HB to be sure we have a full cycle after last HB
lastHB = HB(end);% we keep this to NaN signal after this point
ECG(lastHB:end) = NaN;
HB(end) = [];HB_time(end) = [];

clear HeartBeats
for i = 1:numel(HB)
    HeartBeats(i).R_sample = HB(i);
    HeartBeats(i).R_time = HB_time(i);
end


if hassound
    %% detect heart sounds (S1, S2)
    dataSound = ft_preprocessing(cfg.sound,data);
    HeartSound_bp = dataSound.trial{1};
    soundlag =  round(cfg.sound.tubelength/33000 * data.fsample);% in samples
    HeartSound_bp = [HeartSound_bp(soundlag:end) NaN(1,soundlag-1)];
    
    HeartSoundz = zscore(HeartSound_bp).^2;
    HeartSoundz(lastHB:end) = NaN;
    % SoundPks = peakdetect2(HeartSoundz,1,.1*data.fsample);
    % NB:   peakdetect2 detects peak peak (highest peak above thresh)
    %       peakdetect3 detects onset (first peak above thresh)
    % may want to refine this...
    SoundPks_idx = peakseek(HeartSoundz,.1*data.fsample,1);
    SoundPks_time = time(SoundPks_idx);
    
end

%% Impedance: filtering
dataZ0 = ft_preprocessing(cfg.Z0,data);
Z_0 = - dataZ0.trial{1} * cfg.Z0.factor; % use - by convention
Z_0(lastHB:end) = NaN;

%% Impedance: dZ/dt
dZdt = diff(Z_0);
dZdt = [dZdt NaN]; % pad with one NaN
dZdt = dZdt * data.fsample; % /dt
dZdt(lastHB:end) = NaN;

% [dZdt_max_idx, dZdt_max] = peakdetect2(dZdt,cfg.Z0.factor);
[dZdt_max_idx, dZdt_max] = peakseek(dZdt,0,cfg.Z0.factor);

% [dZdt_min_idx, dZdt_min] = peakdetect2(-dZdt,0,0.1*data.fsample);
[dZdt_min_idx, dZdt_min] = peakseek(-dZdt, 0.11*data.fsample,0);


%% Now compute LVET
% with dZ/dt:
% from time at which dZ/dt reaches cfg.LVE_startfracmax_dZdt * dZ/dt_min (default = 0.15) (negative pointing up)
% to time of following positive peak in dZ/dt

% after each HB, find maximum dZ/dt
alldist = bsxfun(@minus,HB(:),dZdt_max_idx(:)');
dZdt_max_after_HB_idx = NaN(size(HB));dZdt_max_after_HB = NaN(size(HB));dZdt_max_after_HB_time = NaN(size(HB));
for iB = 1:numel(HB)
    try
        j = find(alldist(iB,:) < 0,1,'first');
        dZdt_max_after_HB_idx(iB) = HB(iB) - alldist(iB,j);
        dZdt_max_after_HB_time(iB) = time(dZdt_max_after_HB_idx(iB));
        dZdt_max_after_HB(iB) = dZdt_max(j);
    end
end


[HeartBeats.dZdt_max_value] = rep2struct(dZdt_max_after_HB);
[HeartBeats.dZdt_max_sample] = rep2struct(dZdt_max_after_HB_idx);
[HeartBeats.dZdt_max_time] = rep2struct(dZdt_max_after_HB_time);

% time at which dZdt reaches cfg.LVE_startfracmax_dZdt
dZdt_max_after_HB_15 = dZdt_max_after_HB .* cfg.LVE_startfracmax_dZdt;
LVE_start_time = NaN(size(HB)); LVE_start_dZdt = NaN(size(HB));
for iB = 1:numel(HB)
    x = (HB(iB):dZdt_max_after_HB_idx(iB)) ./ data.fsample;
    dZdt_loc = dZdt(HB(iB):dZdt_max_after_HB_idx(iB));
    i = dsearchn(dZdt_loc(:) - dZdt_max_after_HB_15(iB),0);
    LVE_start_time(iB) = x(i);
    LVE_start_dZdt(iB) = dZdt_loc(i);
end
[HeartBeats.LVE_start_time] = rep2struct(LVE_start_time);
[HeartBeats.LVE_start_sample] = rep2struct(dsearchn(time(:),LVE_start_time(:)));

% after each maximum dZ/dt (after a HB), find minimum dZ/dt
alldist = bsxfun(@minus,dZdt_max_idx(:),dZdt_min_idx(:)');
dZdt_min_after_max_idx = NaN(size(HB));dZdt_min_after_max = NaN(size(HB));dZdt_min_after_max_time = NaN(size(HB));
for iB = 1:numel(HB)
    try
        j = find(alldist(iB,:) < 0,1,'first');
        dZdt_min_after_max_idx(iB) = dZdt_max_after_HB_idx(iB) - alldist(iB,j);
        dZdt_min_after_max(iB) = -dZdt_min(j);
        dZdt_min_after_max_time(iB) = time(dZdt_min_after_max_idx(iB));
    end
end
[HeartBeats.dZdt_min_value] = rep2struct(dZdt_min_after_max);
[HeartBeats.dZdt_min_sample] = rep2struct(dZdt_min_after_max_idx);
[HeartBeats.dZdt_min_time] = rep2struct(dZdt_min_after_max_time);

LVE_stop_time = dZdt_min_after_max_time;
[HeartBeats.LVE_stop_time] = rep2struct(LVE_stop_time);
[HeartBeats.LVE_stop_sample] = rep2struct(dsearchn(time(:),LVE_stop_time(:)));

for iB = 1:numel(HB)
    HeartBeats(iB).Z = mean(Z_0(HeartBeats(iB).LVE_start_sample:HeartBeats(iB).LVE_stop_sample),2);
end

[HeartBeats.LVET] = rep2struct([HeartBeats.LVE_stop_time] - [HeartBeats.LVE_start_time]);



if hassound
    % now with sounds
    % after each HB, pick the first sound
    alldist = bsxfun(@minus,HB(:),SoundPks_idx(:)');
    SoundPks_after_HB_idx = NaN(size(HB));SoundPks_after_HB_time = NaN(size(HB)); SoundPks_after_HB = NaN(size(HB));
    for iB = 1:numel(HB)
        try
            j = find(alldist(iB,:) < 0,1,'first');
            SoundPks_after_HB_idx(iB) = HB(iB) - alldist(iB,j);
            SoundPks_after_HB_time(iB) = time(SoundPks_after_HB_idx(iB));
            SoundPks_after_HB(iB) = HeartSoundz(SoundPks_after_HB_idx(iB));
        end
    end
    
    [HeartBeats.S1_sample] = rep2struct(SoundPks_after_HB_idx);
    [HeartBeats.S1_time] = rep2struct(SoundPks_after_HB_time);
    [HeartBeats.S1_value] = rep2struct(SoundPks_after_HB);
    
    
    alldist = bsxfun(@minus,SoundPks_after_HB_idx(:),SoundPks_idx(:)');
    SoundPks_after_SoundPks_idx = NaN(size(HB));SoundPks_after_SoundPks_time = NaN(size(HB));SoundPks_after_SoundPks = NaN(size(HB));
    for iB = 1:numel(HB)
        try
            j = find(alldist(iB,:) < 0,1,'first');
            SoundPks_after_SoundPks_idx(iB) = SoundPks_after_HB_idx(iB) - alldist(iB,j);
            SoundPks_after_SoundPks_time(iB) = time(SoundPks_after_SoundPks_idx(iB));
            SoundPks_after_SoundPks(iB) = HeartSoundz(SoundPks_after_SoundPks_idx(iB));
        end
    end
    [HeartBeats.S2_sample] = rep2struct(SoundPks_after_SoundPks_idx);
    [HeartBeats.S2_time] = rep2struct(SoundPks_after_SoundPks_time);
    [HeartBeats.S2_value] = rep2struct(SoundPks_after_SoundPks);
end

[HeartBeats.SV] = rep2struct(cfg.SV(cfg.rho,cfg.L,[HeartBeats.Z],[HeartBeats.LVET],[HeartBeats.dZdt_max_value]));
HeartBeats(1).cfg = cfg;

if istrue(cfg.plotallmeas)
    hheart = plot(time,cfg.ECG.magfactor*ECG);
    vline(HB_time,':','color',get(hheart,'color'))
    hold on
    plot(time,Z_0)
    plot(time,dZdt)
    if hassound
        hsound = plot(time,HeartSoundz/10);
        vline(SoundPks_time,':','color',get(hsound,'color'))
    end
    
    xlim([0 5])
    legend({'-Z_0', 'dZ/dt', 'HeartSounds', 'Heart Beats', 'Sound peaks'})
    title('All measures and peaks')
end


if istrue(cfg.interactive)
    hfig = figure(4389309);
    s = get(0,'ScreenSize');
    set(gcf,'position', [1 s(4)/2 s(3) s(4)/2])
    GD = guihandles(hfig);
    GD.time = time;
    GD.ECG = cfg.ECG.magfactor*ECG;
    GD.Z_0 = Z_0;
    GD.dZdt = dZdt;
    GD.HeartSoundz = HeartSoundz;
    GD.HeartBeats = HeartBeats;
    GD.iout = 1;
    GD.LVE_ori = [GD.HeartBeats.LVE_start_time ;GD.HeartBeats.LVE_stop_time]' - [GD.HeartBeats.R_time ;GD.HeartBeats.dZdt_max_time]';
    
    [ioutliers] = find_outliers(GD.LVE_ori);
    [GD.allout,GD.allout_type] = find(ioutliers);
    
    guidata(hfig,GD);
    
    plotoutliers(hfig,[]);
    uiwait(hfig);
end


function plotoutliers(hObject,~)

GD = guidata(hObject);
clf
set(gcf,'Pointer','crosshair')
subplot(1,5,[1 4]);
set(gca,'DefaultLineHittest','off');

hold on

winaround = [-5 5];
switch GD.allout_type(GD.iout)
    case 1
        timearound = winaround + GD.HeartBeats(GD.allout(GD.iout)).LVE_start_time;
    case 2
        timearound = winaround + GD.HeartBeats(GD.allout(GD.iout)).LVE_stop_time;
end
time = clip(GD.time,GD.time,timearound);
HeartBeats = clip(GD.HeartBeats,[GD.HeartBeats.R_time],timearound);

ylim([-30 30])
yl = ylim;

if istrue(GD.HeartBeats(1).cfg.ECG.show)
    ECG = clip(GD.ECG,GD.time,timearound);
    ECG0 =  (yl(2) - diff(yl) / 10);
    hheart = plot(time,ECG + ECG0);
    hline(ECG0,'k');
end

if istrue(GD.HeartBeats(1).cfg.Z0.show)
    Z_0 = clip(GD.Z_0,GD.time,timearound);
    Z_00 = (yl(2) - 2 * diff(yl) / 10);
    plot(time,Z_0 + Z_00);
    hline(Z_00,'k');
end

dZdt = clip(GD.dZdt,GD.time,timearound);
plot(time,dZdt,'tag','dZdt');
vline([HeartBeats.R_time],'color',get(hheart,'color'))

soundfactor = diff(quantile(GD.HeartSoundz,[.001 .999])) / diff(yl)*3;
HeartSoundz = clip(GD.HeartSoundz / soundfactor,GD.time,timearound);
sound0 = yl(1);
hsound = plot(time,HeartSoundz + sound0);
vline([HeartBeats.S1_time],':','color',get(hsound,'color'))
plot([HeartBeats.S1_time],[HeartBeats.S1_value]/ soundfactor + yl(1),'.k')

vline([HeartBeats.S2_time],':','color',get(hsound,'color'))
plot([HeartBeats.S2_time],[HeartBeats.S2_value]/ soundfactor + yl(1),'.k')

vline([HeartBeats.LVE_start_time],'color','k');
plot([HeartBeats.LVE_start_time],GD.dZdt([HeartBeats.LVE_start_sample]),'.k')

vline([HeartBeats.LVE_stop_time],'color','k');
plot([HeartBeats.LVE_stop_time],GD.dZdt([HeartBeats.LVE_stop_sample]),'.k')

switch GD.allout_type(GD.iout)
    case 1
        vline(GD.HeartBeats(GD.allout(GD.iout)).LVE_start_time,'r')
    case 2
        vline(GD.HeartBeats(GD.allout(GD.iout)).LVE_stop_time,'r')
end
xlim(timearound)
xlabel('Time (s)')

set(gca,'ButtonDownFcn',@editpeaks);

subplot(1,5,5)
LVE = [GD.HeartBeats.LVE_start_time ;GD.HeartBeats.LVE_stop_time]' - [GD.HeartBeats.R_time ;GD.HeartBeats.dZdt_max_time]';
hold on;
scatter(LVE(GD.allout(GD.iout),GD.allout_type(GD.iout)),GD.allout_type(GD.iout),'MarkerFaceColor','r');
boxplot(LVE,'plotstyle','compact','orientation','horizontal','jitter',0)
ytick([1 2]);
yticklabel({'HB to max' 'max to min'})


uicontrol('style','pushbutton','string','>',...
    'callback',@next,...
    'units','normalize',...
    'position',[.07 .01 .05 .05])
uicontrol('style','pushbutton','string','<',...
    'callback',@previous,...
    'units','normalize',...
    'position',[.01 .01 .05 .05])
set(gcf,'CloseRequestFcn','tmp = guidata(gcbo);HeartBeats = tmp.HeartBeats;clear tmp;closereq');

function data = clip(data,time,cliptime)
data = data(time>cliptime(1) & time<cliptime(2));


function next(hObject,~)
GD = guidata(hObject);
GD.iout = min(GD.iout+1,numel(GD.allout));
guidata(hObject,GD);
plotoutliers(hObject);

function previous(hObject,~)
GD = guidata(hObject);
GD.iout = max(1,GD.iout-1);
guidata(hObject,GD);
plotoutliers(hObject);


function editpeaks(hObject,~)

GD = guidata(hObject);

pos = get(gca,'CurrentPoint');
p = get(gca,'position');p = p(3:4);
pp = axis;pp = [pp(2)-pp(1) pp(4)-pp(3)];
fac = 1./pp.*1./p;
allpoints = [[GD.HeartBeats.LVE_start_time]' [GD.dZdt([GD.HeartBeats.LVE_start_sample])]';
    [GD.HeartBeats.LVE_stop_time]' [GD.dZdt([GD.HeartBeats.LVE_stop_sample])]'];
HB = repmat(1:numel(GD.HeartBeats),1,2)';

[closest] = dsearchn(bsxfun(@times,allpoints,fac),[pos(1,1:2)].*fac);
iBeat = HB(closest);

dat(1,:) = GD.time;
dat(2,:) = GD.dZdt;
nupoint = dsearchn(bsxfun(@times,dat',fac),pos(1,1:2).*fac);
nudat = dat(:,nupoint);

switch ceil(closest/numel(GD.HeartBeats))
    case 1
        % LVE start
        GD.HeartBeats(iBeat).LVE_start_time = nudat(1);
        GD.HeartBeats(iBeat).LVE_start_sample = nupoint;
    case 2
        % LVE end
        GD.HeartBeats(iBeat).LVE_stop_time = nudat(1);
        GD.HeartBeats(iBeat).LVE_stop_sample = nupoint;
end

guidata(hObject,GD);
plotoutliers(hObject);
