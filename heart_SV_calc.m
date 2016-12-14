function HeartBeats = heart_SV_calc(cfg,data,HeartBeats)

% [HeartBeats] = heart_SV_calc(cfg)
% [HeartBeats] = heart_SV_calc(cfg,data)
% [HeartBeats] = heart_SV_calc(cfg,data,HeartBeats)
%
%   The results are returned in HeartBeats, a vector structure with as many
%   elements as detected R peaks in the ECG. Retrieve (e.g.) Stroke volume
%   for each heart beat with : SV = [HeartBeats.SV] (square brackets are
%   important)
% 
% If you are calling heart_SV_calc with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%         cfg.dataset       = string with the filename
%
% If you are calling heart_SV_calc with also the second input argument
% "data", then that should contain data that was read from file in a
% previous call to FT_PREPROCESSING. The configuration options below apply.
%
% If HeartBeats is provided as a third argument, it is assumed to be the
% output of heart_peak_detect. No detection is computed on the ECG. If a SV
% field is present (as produced by heart_SV_calc, this function) in
% HeartBeats, then only an interactive plot is opened.
%
% The options are specified with:
%     - ECG options:
%         cfg.ECG.channel     = channel to use for ECG (default = 'ECG1')
%         cfg.ECG.show        = whether to show the ECG on any plot (default = 'yes')
%         cfg.ECG_peakcfg     = options to heart_peak_detect (see help heart_peak_detect) (default = [], use defaults)
%     - Impedance (Z0) channel options
%         cfg.Z0.channel      = channel name to use for Z0 (default = 'EBI_-_Magnitude')
%         cfg.Z0.lpfilter     = whether or not to low pass filter Z0 data (default = 'yes')
%         cfg.Z0.lpfreq       = frequency cutoff for low pass filter (default = 15)
%         cfg.Z0.lpfilttype   = low pass filter type (default = 'but')
%         cfg.Z0.factor       = impedance factor in ohms/volt (default = 5)
%         cfg.Z0.show         = whether or not to show Z0 in plots (default = 'yes')
%     - Sound channel options:
%         cfg.sound.channel   = channel name to use (default = 'HeartSound')
%         cfg.sound.bpfilter  = whether to band pass filter the sound channel (default = 'yes')
%         cfg.sound.bpfreq    = pass band frequency in Hz (default = [40 60])
%         cfg.sound.tubelength = distance in cm the sound has to travel before reaching the microphone (default = 152, --> delay ~ 4 ms)
% 
%     - Algorithmic limits and Stroke volume computation constants:
%         cfg.QCmax           = maximal duration of the QC interval in seconds (default = .3)
%         cfg.RS1max          = maximal time interval between R peak and S1 in seconds (default = .2)
%         cfg.LVETmax         = maximal left ventricular ejection time in seconds (default = .4)
%         cfg.LVE_startfracmax_dZdt = fraction of maximum dZdt value to use to detect ejection start (default = 0.15, Kubicek, W. G., R. P. Patterson, and D. A. Witsoe. “Impedance Cardiography as a Noninvasive Method of Monitoring Cardiac Function and Other Parameters of the Cardiovascular System*.” Annals of the New York Academy of Sciences 170, no. 2 (July 1, 1970): 724–32)
%         cfg.rho             = resistivity of blood in Ohms (default = 135, see Berntson, Gary G., Karen S. Quigley, and Dave Lozano. “Cardiovascular Psychophysiology.” Handbook of Psychophysiology 3 (2007): 182–210.)
%         cfg.L               = distance between the two middle electrodes in cm (default = 30)
% 
%     - Ploting options:
%         cfg.interactive     = whether or not an interactive plot should be opened at the end to allow changint outlier peaks (default = 'yes')
%         cfg.ECG.scale       = display scale for ECG   (default = should look ok)
%         cfg.Z0.scale        = display scale for Z0    (default = should look ok)
%         cfg.dZdt.scale      = display scale for dZdt  (default = should look ok)
%         cfg.sound.scale     = display scale for sound (default = should look ok)
%
%
%   Algorithm:
%
%         Peaks are first detected on the ECG signal using heart_peak_detect with
%         the options specified. The last heart beat is discarded to ensure all
%         heart cycles are complete. The impedance channel is low pass filtered,
%         dZ/dt is computed.
%         For each detected heart cycle, start with the Q point. Find the maximum
%         dZ/dt (dZdt_max) following this Q point (within QCmax seconds). Find the time at
%         which dZ/dt reaches LVE_startfracmax_dZdt * dZdt_max. This is the B
%         point (start of the left ventricular ejection (LVE)). Find the time
%         following dZdt_max where dZ/dt reaches a minimum value (within LVETmax
%         seconds of the B point). This is the X point (end of LVE).
%         Compute LVET, the duration in seconds of LVE.
%         If heart sounds are provided, the sound signal is band pass filtered
%         according to cfg.sound.bpfilter, zscored and squared. S1 and S2 are
%         detected as the maxima of this transformed signal. S1 is the max between
%         R peak and RS1max, S2 is the max between S1 and LVETmax. After S1 and S2
%         have been determined, a check for outliers on X-S2 duration is performed.
%         For each outlier, if a negative peak in dZdt is present between the
%         previously found X point (negative peak after dZdt_max) and S2, this new
%         peak is taken as X point. Otherwise leave as is.
%         Compute Z, the mean Z0 during LVE.
%         Compute Stroke volume with the following formula:
% 
%               SV = rho .* L .^ 2 ./ Z .^ 2 .* LVET .* dZdt_max
% 

% v0 Maximilien Chaumon November 2016


def             = [];
def.ECG.channel = 'ECG1';
def.ECG.scale   = [];
def.ECG.show    = 'yes';

def.ECG_peakcfg = [];
def.ECG_peakcfg.channel = def.ECG.channel;

def.Z0.channel  = 'EBI_-_Magnitude';
def.Z0.scale    = [];
def.Z0.lpfilter = 'yes';
def.Z0.lpfreq   = 15;
def.Z0.lpfilttype = 'but';
def.Z0.show     = 'yes';

def.dZdt.scale  = [];

def.sound.channel = 'HeartSound';
def.sound.scale   = [];
def.sound.bpfilter = 'yes';
def.sound.bpfreq = [40 60];
def.sound.tubelength = 152;

def.QCmax       = .3;
def.RS1max      = .2;
def.LVETmax     = .4;
def.LVE_startfracmax_dZdt = 0.15;
def.rho         = 135;
def.L           = 30;
def.SV          = @(rho,L,Z,LVET,dZdt_max) rho .* L .^ 2 ./ Z .^ 2 .* LVET .* dZdt_max;

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

cfg.ECG_peakcfg.fsample = data.fsample;

if ~exist('HeartBeats','var')
    [HeartBeats] = heart_peak_detect(cfg.ECG_peakcfg,data);
end
% remove last HB to be sure we have a full cycle after last HB
lastHB = HeartBeats(end);% we keep this to NaN signal after this point
ECG(lastHB.R_sample:end) = NaN;
HeartBeats(end) = [];
HB = [HeartBeats.R_sample];
if isempty(cfg.ECG.scale)
    cfg.ECG.scale = 1/(diff(quantile(ECG(:),[.001 .999])) / 10) * sign(skewness(ECG));
end


if hassound
    %% detect heart sounds (S1, S2)
    dataSound = ft_preprocessing(cfg.sound,data);
    HeartSound_bp = dataSound.trial{1};
    soundlag =  round(cfg.sound.tubelength/33000 * data.fsample);% in samples
    HeartSound_bp = [HeartSound_bp(soundlag:end) NaN(1,soundlag-1)];

    HeartSoundz = nanzscore(HeartSound_bp).^2;
    HeartSoundz(lastHB.R_sample:end) = NaN;
    if isempty(cfg.sound.scale)
        cfg.sound.scale = 1/(diff(quantile(HeartSoundz(:),[.001 .999])) / 10);
    end
    % SoundPks = peakdetect2(HeartSoundz,1,.1*data.fsample);
    % NB:   peakdetect2 detects peak peak (highest peak above thresh)
    %       peakdetect3 detects onset (first peak above thresh)
    % may want to refine this...
    SoundPks_idx = peakseek(HeartSoundz,1,.1*data.fsample);
    SoundPks_time = time(SoundPks_idx);
else
    HeartSoundz = [];
end

%% Impedance: filtering
dataZ0 = ft_preprocessing(cfg.Z0,data);
Z0 = - dataZ0.trial{1}; % use - by convention
Z0(lastHB.R_sample:end) = NaN;
if isempty(cfg.Z0.scale)
    cfg.Z0.scale = 1/(diff(quantile(Z0(:),[.001 .999])) / 10);
end


%% Impedance: dZ/dt
dZdt = diff(Z0);
dZdt = [dZdt NaN]; % pad with one NaN
dZdt = dZdt * data.fsample; % /dt
dZdt(lastHB.R_sample:end) = NaN;
if isempty(cfg.dZdt.scale)
    cfg.dZdt.scale = 1/(diff(quantile(Z0(:),[.001 .999])) / 3);
end
if isfield(HeartBeats,'SV')
    % consider we have run this already. Just plot and leave
    uiwait(plotsetup(cfg,time,ECG,Z0,dZdt,HeartSoundz,HeartBeats))
    tmp = guidata(gcf);HeartBeats = tmp.HeartBeats;
    delete(gcf)
    return
end

%% now find points and compute LVET
for i_R = 1:numel(HeartBeats)
    % start from Q peak
    Q = HeartBeats(i_R).Q_sample;
    idx = Q:Q + cfg.QCmax * data.fsample;
    dZdttmp = dZdt(idx);
    % find dZdtmax
    [v,p] = max(dZdttmp);
    dZdt_max_value(i_R) = v;
    dZdt_max_sample(i_R) = p + idx(1) - 1;
    dZdt_max_time(i_R) = time(dZdt_max_sample(i_R));
    % deduce B point (PEP)
    idx = Q:dZdt_max_sample(i_R);
    dZdttmp = dZdt(idx);
    dZdt_end_PEP = dZdt_max_value(i_R) .* cfg.LVE_startfracmax_dZdt;
    [v,p] = min(abs(dZdttmp - dZdt_end_PEP));
    B_sample(i_R) = p + idx(1) - 1;
    B_time(i_R) = time(B_sample(i_R));
    B_value(i_R) = v;
    PEP_sample(i_R) = B_sample(i_R) - Q;
    PEP_time(i_R) = PEP_sample(i_R) / data.fsample;
    % find dZdtmin (X point)
    idx = B_sample(i_R):B_sample(i_R) + cfg.LVETmax * data.fsample;
    dZdttmp = dZdt(idx);
    [v,p] = min(dZdttmp);
    X_value(i_R) = v;
    X_sample(i_R) = idx(1) + p - 1;
    X_time(i_R) = time(X_sample(i_R));
    % deduce LVET
    LVET_sample(i_R) = X_sample(i_R) - B_sample(i_R);
    LVET_time(i_R) = X_time(i_R) - B_time(i_R);
    Z(i_R) = mean(Z0(B_sample(i_R):X_sample(i_R)),2);
    
end
[HeartBeats.dZdt_max_value] = rep2struct(dZdt_max_value);
[HeartBeats.dZdt_max_sample] = rep2struct(dZdt_max_sample);
[HeartBeats.dZdt_max_time] = rep2struct(dZdt_max_time);
[HeartBeats.B_sample] = rep2struct(B_sample);
[HeartBeats.B_time] = rep2struct(B_time);
[HeartBeats.B_value] = rep2struct(B_value);
[HeartBeats.PEP_sample] = rep2struct(PEP_sample);
[HeartBeats.PEP_time] = rep2struct(PEP_time);
[HeartBeats.X_sample] = rep2struct(X_sample);
[HeartBeats.X_time] = rep2struct(X_time);
[HeartBeats.X_value] = rep2struct(X_value);
[HeartBeats.LVET_sample] = rep2struct(LVET_sample);
[HeartBeats.LVET_time] = rep2struct(LVET_time);
[HeartBeats.Z] = rep2struct(Z);

if hassound
    % now with sounds
    % after each HB, pick the first sound
    for i_R = 1:numel(HeartBeats)
        idx = HeartBeats(i_R).R_sample:HeartBeats(i_R).R_sample + cfg.RS1max * data.fsample;
        tmp = HeartSoundz(idx);
        [v,p] = max(tmp);
        S1_sample(i_R) = idx(1) + p - 1;
        S1_time(i_R) = time(S1_sample(i_R));
        S1_value(i_R) = v;
    end
    [HeartBeats.S1_sample] = rep2struct(S1_sample);
    [HeartBeats.S1_time] = rep2struct(S1_time);
    [HeartBeats.S1_value] = rep2struct(S1_value);

    % after each S1, pick S2
    for i_R = 1:numel(HeartBeats)
        idx = HeartBeats(i_R).S1_sample:HeartBeats(i_R).S1_sample + cfg.LVETmax * data.fsample;
        tmp = HeartSoundz(idx);
        [v,p] = max(tmp);
        S2_sample(i_R) = idx(1) + p - 1;
        S2_time(i_R) = time(S2_sample(i_R));
        S2_value(i_R) = v;
    end
    [HeartBeats.S2_sample] = rep2struct(S2_sample);
    [HeartBeats.S2_time] = rep2struct(S2_time);
    [HeartBeats.S2_value] = rep2struct(S2_value);
    
    % try an automatic fix of outliers based on bad X point
    outliers = find(find_outliers([HeartBeats.S2_time]-[HeartBeats.X_time]));
    for i_out = 1:numel(outliers)
        idx = HeartBeats(outliers(i_out)).X_sample:HeartBeats(outliers(i_out)).S2_sample;
        tmp = dZdt(idx);
        % if there is one or more peaks in dZdt before S2, we take the last
        % one as our new X point.
        [p,v] = peakseek(-tmp,-Inf,0);
        
        if not(isempty(p))
            HeartBeats(outliers(i_out)).X_sample = p(end) + idx(1) - 1;
            HeartBeats(outliers(i_out)).X_time = time(HeartBeats(outliers(i_out)).X_sample);
            HeartBeats(outliers(i_out)).X_value = -v(end);
            HeartBeats(outliers(i_out)).LVET_sample = HeartBeats(outliers(i_out)).X_sample - HeartBeats(outliers(i_out)).B_sample;
            HeartBeats(outliers(i_out)).LVET_time = HeartBeats(outliers(i_out)).X_time - HeartBeats(outliers(i_out)).B_time;
            HeartBeats(outliers(i_out)).Z = mean(Z0(HeartBeats(outliers(i_out)).B_sample:HeartBeats(outliers(i_out)).X_sample),2);
        end
    end
    
end

[HeartBeats.SV] = rep2struct(cfg.SV(cfg.rho,cfg.L,[HeartBeats.Z],[HeartBeats.LVET_time],[HeartBeats.dZdt_max_value]));
HeartBeats(1).cfg = cfg;

if istrue(cfg.interactive)
    uiwait(plotsetup(cfg,time,ECG,Z0,dZdt,HeartSoundz,HeartBeats))
    tmp = guidata(gcf);HeartBeats = tmp.HeartBeats;
    delete(gcf)
    return
end

function hfig = plotsetup(cfg,time,ECG,Z0,dZdt,HeartSoundz,HeartBeats)

hfig = figure(4389309);

s = get(0,'ScreenSize');
set(gcf,'position', [1 s(4)/2 s(3) s(4)/2])
GD = guihandles(hfig);
GD.axtimecourse = subplot(1,5,[1 4]);
GD.axB = subplot(3,5,5);
GD.axX = subplot(3,5,10);
GD.axLVET = subplot(3,5,15);
GD.time = time;
GD.winaround = [-2.5 2.5];
GD.ECG = ECG;
GD.ECGscale = cfg.ECG.scale;
GD.Z0 = Z0;
GD.Z0scale = cfg.Z0.scale;
GD.dZdt = dZdt;
GD.dZdtscale = cfg.dZdt.scale;
if ~isempty(HeartSoundz)
    GD.HeartSoundz = HeartSoundz;
    GD.HeartSoundzscale = cfg.sound.scale;
end
GD.HeartBeats = HeartBeats;
GD.iout = 1;
GD.LVE_ori = [GD.HeartBeats.B_time ;GD.HeartBeats.X_time]' - [GD.HeartBeats.R_time ;GD.HeartBeats.dZdt_max_time]';

[ioutliers] = find_outliers(GD.LVE_ori);
[GD.allout,GD.allout_type] = find(ioutliers);

guidata(hfig,GD);
set(gcf,'CloseRequestFcn','uiresume');
set(gcf,'Pointer','crosshair')

plotoutliers(hfig,[]);


function plotoutliers(hObject,~)

% retrieve data
GD = guidata(hObject);

% work on time course axis
axes(GD.axtimecourse);
cla
set(gca,'DefaultLineHittest','off');
cols = get(gca,'colororder');

hold on
ytl = {};
yt = [];

% set up zoom to time window around current outlier
switch GD.allout_type(GD.iout)
    case 1
        timearound = GD.winaround + GD.HeartBeats(GD.allout(GD.iout)).B_time;
    case 2
        timearound = GD.winaround + GD.HeartBeats(GD.allout(GD.iout)).X_time;
end

% clip data to time window of interest
time = clip(GD.time,GD.time,timearound);
HeartBeats = clip(GD.HeartBeats,[GD.HeartBeats.R_time],timearound);

ylim([-30 30])
yl = ylim;

% show ECG
if istrue(GD.HeartBeats(1).cfg.ECG.show)
    ECG = clip(GD.ECG,GD.time,timearound);
    ECG0 =  (yl(2) - diff(yl) / 10);
    hheart = plot(time,ECG * GD.ECGscale + ECG0,'color',cols(1,:));
    hline(ECG0,':k');
    yt(end+1) = ECG0;
    ytl{end+1} = 'ECG';
end
% show Z0
if istrue(GD.HeartBeats(1).cfg.Z0.show)
    Z0 = clip(GD.Z0,GD.time,timearound);
    Z0 = Z0 - mean(Z0);
    Z00 = (yl(2) - 2 * diff(yl) / 10);
    plot(time,Z0 * GD.Z0scale + Z00,'color',cols(2,:));
    hline(Z00,':k');
    yt(end+1) = Z00;
    ytl{end+1} = 'Z_0';
end
% show dZ/dt
dZdt = clip(GD.dZdt,GD.time,timearound);
plot(time,dZdt * GD.dZdtscale,'tag','dZdt','color',cols(3,:));
hline(0,':k')
yt(end+1) = 0;
ytl{end+1} = 'dZ/dt';

% add all R peaks same color as ECG
vline([HeartBeats.R_time],'color',get(hheart,'color'))

% add sound with peaks
if isfield(GD,'HeartSoundz')
    HeartSoundz = clip(GD.HeartSoundz,GD.time,timearound);
    sound0 = yl(1);
    hsound = plot(time,HeartSoundz * GD.HeartSoundzscale + sound0,'color',cols(4,:));
    yt(end+1) = sound0;
    ytl{end+1} = 'Heart Sound';
    vline([HeartBeats.S1_time],':','color',get(hsound,'color'))
    plot([HeartBeats.S1_time],[HeartBeats.S1_value] * GD.HeartSoundzscale + sound0,'.k')

    vline([HeartBeats.S2_time],':','color',get(hsound,'color'))
    plot([HeartBeats.S2_time],[HeartBeats.S2_value] * GD.HeartSoundzscale + sound0,'.k')
end

% add B points
vline([HeartBeats.B_time],'color','k');
plot([HeartBeats.B_time],GD.dZdt([HeartBeats.B_sample]) * GD.dZdtscale,'.k')

% add X points
vline([HeartBeats.X_time],'color','k');
plot([HeartBeats.X_time],GD.dZdt([HeartBeats.X_sample]) * GD.dZdtscale,'.k')

% highlight current outlier
switch GD.allout_type(GD.iout)
    case 1
        vline(GD.HeartBeats(GD.allout(GD.iout)).B_time,'r')
    case 2
        vline(GD.HeartBeats(GD.allout(GD.iout)).X_time,'r')
end

xlim(timearound)
xlabel('Time (s)')
ytick(yt(end:-1:1))
yticklabel(ytl(end:-1:1))

set(gca,'ButtonDownFcn',@editpeaks);
title(sprintf('Outlier %d / %d',GD.iout,numel(GD.allout)));

% add boxplots
ax = [GD.axB GD.axX];
LVE = [GD.HeartBeats.B_time ;GD.HeartBeats.X_time]' - [GD.HeartBeats.R_time ;GD.HeartBeats.dZdt_max_time]';

% boxplot for B points
axes(GD.axB)
if isempty(findobj(gca,'tag','Box'))
    boxplot(LVE(:,1),'plotstyle','compact','orientation','horizontal','symbol','')
end
hold on;
ytick([1]);
title('R peak to B point')

% boxplot for X points
axes(GD.axX)
if isempty(findobj(gca,'tag','Box'))
    boxplot(LVE(:,2),'plotstyle','compact','orientation','horizontal','symbol','')
end
hold on;
title('dZdt_{max} to X point')

% add all outliers
try delete(GD.scatout); end
try delete(GD.scatiout); end
    
axes(ax(GD.allout_type(GD.iout)));
GD.scatout = scatter(LVE(GD.allout(GD.iout),GD.allout_type(GD.iout)),1,'MarkerFaceColor','r','MarkerEdgeColor','none');
for iout = 1:numel(GD.allout)
    cb = {@go,iout};
    axes(ax(GD.allout_type(iout)));
    GD.scatiout(iout) = scatter(LVE(GD.allout(iout),GD.allout_type(iout)),1,'MarkerEdgeColor','b','ButtonDownFcn',cb);
end

% add LVET boxplot
LVET = diff(LVE,[],2);
axes(GD.axLVET);
cla
boxplot(LVET,'plotstyle','compact','orientation','horizontal','jitter',0)
hold on;
scatter(LVET(GD.allout(GD.iout)),1,'MarkerFaceColor','r','MarkerEdgeColor','none');
title('LVET')
xlabel('Time (s)')

% store data
guidata(gcf,GD);


% add buttons
uicontrol('style','pushbutton','string','<',...
    'callback',@previous,...
    'units','normalize',...
    'position',[.01 .01 .05 .05])
uicontrol('style','pushbutton','string','>',...
    'callback',@next,...
    'units','normalize',...
    'position',[.07 .01 .05 .05])
uicontrol('style','pushbutton','string','zoom -',...
    'callback',@zminus,...
    'units','normalize',...
    'position',[.13 .01 .05 .05])
uicontrol('style','pushbutton','string','zoom +',...
    'callback',@zplus,...
    'units','normalize',...
    'position',[.19 .01 .05 .05])

function data = clip(data,time,cliptime)
% a simple function to clip a specific time window of data
data = data(time>cliptime(1) & time<cliptime(2));

function zplus(hObject,~)
% zoom in
GD = guidata(hObject);
GD.winaround = GD.winaround / 2;
guidata(hObject,GD);
plotoutliers(hObject);

function zminus(hObject,~)
% zoom out
GD = guidata(hObject);
GD.winaround = GD.winaround * 2;
guidata(hObject,GD);
plotoutliers(hObject);

function go(hObject,~,iout)
% step to specific outlier (click in boxplots)
GD = guidata(hObject);
GD.iout = iout;
guidata(hObject,GD);
plotoutliers(hObject);

function next(hObject,~)
% step to next outlier
GD = guidata(hObject);
GD.iout = min(GD.iout+1,numel(GD.allout));
guidata(hObject,GD);
plotoutliers(hObject);

function previous(hObject,~)
% step to previous outlier
GD = guidata(hObject);
GD.iout = max(1,GD.iout-1);
guidata(hObject,GD);
plotoutliers(hObject);


function editpeaks(hObject,~)
% click on time course figure

% retrieve data
GD = guidata(hObject);

% where did we click?
pos = get(gca,'CurrentPoint');
% try to equalize distances along x and y axes
p = get(gca,'position');p = p(3:4);
pp = axis;pp = [pp(2)-pp(1) pp(4)-pp(3)];
fac = 1./pp.*1./p;

% current beat
iBeat = GD.allout(GD.iout);

dat(1,:) = GD.time;
dat(2,:) = GD.dZdt;
nupoint = dsearchn(bsxfun(@times,dat',fac),pos(1,1:2).*fac);
nudat = dat(:,nupoint);

switch GD.allout_type(GD.iout)%ceil(closest/numel(GD.HeartBeats))
    case 1
        % LVE start
        GD.HeartBeats(iBeat).B_time = nudat(1);
        GD.HeartBeats(iBeat).B_sample = nupoint;
    case 2
        % LVE end
        GD.HeartBeats(iBeat).X_time = nudat(1);
        GD.HeartBeats(iBeat).X_sample = nupoint;
end
GD.HeartBeats(iBeat).LVET_sample = GD.HeartBeats(iBeat).X_sample - GD.HeartBeats(iBeat).B_sample;
GD.HeartBeats(iBeat).LVET_time = GD.HeartBeats(iBeat).X_time - GD.HeartBeats(iBeat).B_time;
GD.HeartBeats(iBeat).Z = mean(GD.Z0(GD.HeartBeats(iBeat).B_sample:GD.HeartBeats(iBeat).X_sample),2);
cfg = GD.HeartBeats(1).cfg;
GD.HeartBeats(iBeat).SV = cfg.SV(cfg.rho,cfg.L,[GD.HeartBeats(iBeat).Z],[GD.HeartBeats(iBeat).LVET_time],[GD.HeartBeats(iBeat).dZdt_max_value]);

guidata(hObject,GD);
plotoutliers(hObject);

function [z,mu,sigma] = nanzscore(x,flag,dim)

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = x; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

z = NaN(size(x));
inan = ~isnan(x);
x = x(inan);
[ztmp,mu,sigma] = zscore(x,flag,dim);
z(inan) = ztmp;


function [locs, pks]=peakseek(x,minpeakh,minpeakdist)
% [locs, pks]=peakseek(x,minpeakh,minpeakdist)
% Alternative to the findpeaks function.  This thing runs much much faster.
% It really leaves findpeaks in the dust.  It also can handle ties between
% peaks.  Findpeaks just erases both in a tie.  Shame on findpeaks.
%
% x is a vector input (generally a timecourse)
% minpeakdist is the minimum desired distance between peaks (optional, defaults to 1)
% minpeakh is the minimum height of a peak (optional)
%
% (c) 2010
% Peter O'Connor
% peter<dot>ed<dot>oconnor .AT. gmail<dot>com

if size(x,2)==1, x=x'; end

% Find all maxima and ties
locs=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1;

if nargin<2, minpeakdist=1; end % If no minpeakdist specified, default to 1.

if nargin>2 % If there's a minpeakheight
    locs(x(locs)<=minpeakh)=[];
end

if minpeakdist>1
    while 1

        del=diff(locs)<minpeakdist;

        if ~any(del), break; end

        pks=x(locs);

        [garb mins]=min([pks(del) ; pks([false del])]); %#ok<ASGLU>

        deln=find(del);

        deln=[deln(mins==1) deln(mins==2)+1];

        locs(deln)=[];

    end
end

if nargout>1,
    pks=x(locs);
end
function sk = skewness(data)
x0 = data - nanmean(data);
s2 = nanmean(data.^2);
m3 = nanmean(data.^3);
sk = m3 ./ s2 .^(1.5);
function [nb,channame,strnames] = chnb(channame, varargin)

% chnb() - return channel number corresponding to channel names in an EEG
%           structure
%
% Usage:
%   >> [nb]                 = chnb(channameornb);
%   >> [nb,names]           = chnb(channameornb,...);
%   >> [nb,names,strnames]  = chnb(channameornb,...);
%   >> [nb]                 = chnb(channameornb, EEG);
%   >> [nb]                 = chnb(channameornb, labels);
%
% Input:
%   channameornb  - If a string or cell array of strings, it is assumed to
%                   be (part of) the name of channels to search. Either a
%                   string with space separated channel names, or a cell
%                   array of strings. 
%                   Note that regular expressions can be used to match
%                   several channels. See regexp.
%                   If only one channame pattern is given and the string
%                   'inv' is attached to it, the channels NOT matching the
%                   pattern are returned.
%   labels        - Channel names as found in {EEG.chanlocs.labels}.
%
% Output:
%   nb            - Channel numbers in labels, or in the EEG structure
%                   found in the caller workspace (i.e. where the function
%                   is called from) or in the base workspace, if no EEG
%                   structure exists in the caller workspace.
%   names         - Channel names, cell array of strings.
%   strnames      - Channel names, one line character array.
narginchk(1,2);
if nargin == 2
    if isstruct(varargin{1}) && isfield(varargin{1},'setname')
        % assume it's an EEG dataset
        labels = {varargin{1}(1).chanlocs.labels};
    else
        labels = varargin{1};
    end
else
    
    try
        EEG = evalin('caller','EEG');
    catch
        try
            EEG = evalin('base','EEG');
        catch
            error('Could not find EEG structure');
        end
    end
    if not(isfield(EEG,'chanlocs'))
        error('No channel list found');
    end
    EEG = EEG(1);
    labels = {EEG.chanlocs.labels};
end
if iscell(channame) || ischar(channame)
    
    if ischar(channame) || iscellstr(channame)
        if iscellstr(channame) && numel(channame) == 1 && isempty(channame{1})
            channame = '';
        end
        tmp = regexp(channame,'(\S*) ?','tokens');
        channame = {};
        for i = 1:numel(tmp)
            if iscellstr(tmp{i}{1})
                channame{i} = tmp{i}{1}{1};
            else
                channame{i} = tmp{i}{1};
            end
        end
        if isempty(channame)
            nb = [];
            channame = {};
            strnames = '';
            return
        end
    end
    if numel(channame) == 1 && not(isempty(strmatch('inv',channame{1})))
        cmd = 'exactinv';
        channame{1} = strrep(channame{1},'inv','');
    else
        channame{1} = channame{1};
        cmd = 'exact';
    end
    nb = regexpcell(labels,channame,[cmd 'ignorecase']);
    
elseif isnumeric(channame)
    nb = channame;
    if nb > numel(labels)
        nb = [];
    end
end
channame = labels(nb);
strnames = sprintf('%s ',channame{:});
if not(isempty(strnames))
    strnames(end) = [];
end
if nargout == 0
    disp(channame)
    disp(nb)
    clear
end


function idx = regexpcell(c,pat, cmds)

% idx = regexpcell(c,pat, cmds)
%
% Return indices idx of cells in c that match pattern(s) pat (regular expression).
% Pattern pat can be char or cellstr. In the later case regexpcell returns
% indexes of cells that match any pattern in pat.
%
% cmds is a string that can contain one or several of these commands:
% 'inv' return indexes that do not match the pattern.
% 'ignorecase' will use regexpi instead of regexp
% 'exact' performs an exact match (regular expression should match the whole strings in c).
% 'all' (default) returns all indices, including repeats (if several pat match a single cell in c).
% 'unique' will return unique sorted indices.
% 'intersect' will return only indices in c that match ALL the patterns in pat.
% 
% v1 Maximilien Chaumon 01/05/09
% v1.1 Maximilien Chaumon 24/05/09 - added ignorecase
% v2 Maximilien Chaumon 02/03/2010 changed input method.
%       inv,ignorecase,exact,combine are replaced by cmds

narginchk(2,3)
if not(iscellstr(c))
    error('input c must be a cell array of strings');
end
if nargin == 2
    cmds = '';
end
if not(isempty(regexpi(cmds,'inv', 'once' )))
    inv = true;
else
    inv = false;
end
if not(isempty(regexpi(cmds,'ignorecase', 'once' )))
    ignorecase = true;
else
    ignorecase = false;
end
if not(isempty(regexpi(cmds,'exact', 'once' )))
    exact = true;
else
    exact = false;
end
if not(isempty(regexpi(cmds,'unique', 'once' )))
    combine = 2;
elseif not(isempty(regexpi(cmds,'intersect', 'once' )))
    combine = 3;
else
    combine = 1;
end

if ischar(pat)
    pat = cellstr(pat);
end

if exact
    for i_pat = 1:numel(pat)
        pat{i_pat} = ['^' pat{i_pat} '$'];
    end
end
    
for i_pat = 1:length(pat)
    if ignorecase
        trouv = regexpi(c,pat{i_pat}); % apply regexp on each pattern
    else
        trouv = regexp(c,pat{i_pat}); % apply regexp on each pattern
    end
    idx{i_pat} = [];
    for i = 1:numel(trouv)
        if not(isempty(trouv{i}))% if there is a match, store index
            idx{i_pat}(end+1) = i;
        end
    end
end
switch combine
    case 1
        idx = [idx{:}];
    case 2
        idx = unique([idx{:}]);
    case 3
        for i_pat = 2:length(pat)
            idx{1} = intersect(idx{1},idx{i_pat});
        end
        idx = idx{1};
end
if inv % if we want to invert result, then do so.
    others = 1:numel(trouv);
    others(idx) = [];
    idx = others;
end


function [i_outliers]= find_outliers(data,w)

% [i_outliers]= find_outliers(data,w)
% returns logical yes for outliers in data.
% outliers are defined as data points beyond [q(1) - w * IQR, q(2) + w * IQR]
% with 
%   [q] = quantile(data,[0.25 0.75])
%   IQR = diff(q)
% default value for w is 1.5
% 
% if data is a matrix, separate columns are processed independently.
%
% Max 2016

narginchk(1,2);
if nargin == 1
    w = 1.5;
end
if isvector(data)
    data = data(:);
end
i_outliers = false(size(data));
for i = 1:size(data,2)
    q = quantile(data(:,i),[.25 .75]);
    outlim = q + [-w w] .* diff(q);
    i_outliers(:,i) = data(:,i) < outlim(1) |data(:,i) > outlim(2);
end
function h = vline(x,varargin)

% h = vline(x,varargin)
% add vertical line(s) on the current axes at x
% optional inputs: 
%   'ticks', [numeric] : position of ticks on the line
%   'ticklength', [scalar]: tick length nth of size of the axis
%                           default = 1/40: ticklength = diff(ylim)/40 
%   'color', [numeric] : color of the lines. may have as many rows as there are lines.
% all other varargin arguments are passed to plot...

x = x(:);
varg = cellfun(@(x)num2str(x),varargin,'uniformoutput',0);
ticks = cellfun(@(x)strcmp(x,'ticks'),varg);
ticklength = cellfun(@(x)strcmp(x,'ticklength'),varg);
if any(ticks)
    iticks = find(ticks);
    ticks = varargin{iticks+1};
    if any(ticklength)
        iticklength = find(ticklength);
        ticklength = varargin{iticklength+1};
        varargin(iticklength:iticklength+1) = [];
    else
        ticklength = 1/40;
    end
    varargin(iticks:iticks+1) = [];
    xs = [x - diff(xlim)*ticklength, x + diff(xlim)*ticklength];
    arrayfun(@(i) line(xs,repmat(ticks(i),size(xs,1),2),'color','k'),1:numel(ticks));
end
ho = ishold;
hold on
c = cellfun(@(x)strcmp(x,'color'),varg);
if any(c)
    cs = varargin{find(c)+1};
    varargin([find(c),find(c)+1]) = [];
end
h = plot([x x]',repmat(ylim,numel(x),1)',varargin{:});
if any(c)
    if numel(h) == size(cs,1)
        for ih = 1:numel(h)
            set(h(ih),'color',cs(ih,:))
        end
    else
        set(h,'color',cs(1,:))
    end
end
if not(ho)
    hold off
end
if nargout == 0
    clear h
end
function h = hline(y,varargin)

% h = hline(y,varargin)
% add horizontal line(s) on the current axes at y
% optional inputs: 
%   'ticks', [numeric] : position of ticks on the line
%   'ticklength', [scalar]: tick length nth of size of the axis
%                           default = 1/40: ticklength = diff(ylim)/40 
%   'color', [numeric] : color of the lines. may have as many rows as there are lines.
% all other varargin arguments are passed to plot...

y = y(:);
varg = cellfun(@(x)num2str(x),varargin,'uniformoutput',0);
ticks = cellfun(@(x)strcmp(x,'ticks'),varg);
ticklength = cellfun(@(x)strcmp(x,'ticklength'),varg);
if any(ticks)
    iticks = find(ticks);
    ticks = varargin{iticks+1};
    if any(ticklength)
        iticklength = find(ticklength);
        ticklength = varargin{iticklength+1};
        varargin(iticklength:iticklength+1) = [];
    else
        ticklength = 1/40;
    end
    varargin(iticks:iticks+1) = [];
    ys = [y - diff(ylim)*ticklength, y + diff(ylim)*ticklength];
    arrayfun(@(i) line(repmat(ticks(i),size(ys,1),2),ys,'color','k'),1:numel(ticks));
end
ho = ishold;
hold on
c = cellfun(@(x)strcmp(x,'color'),varg);
if any(c)
    cs = varargin{find(c)+1};
    varargin([find(c),find(c)+1]) = [];
end
h = plot(repmat(xlim,numel(y),1)',[y y]',varargin{:});
if any(c)
    if numel(h) == size(cs,1)
        for ih = 1:numel(h)
            set(h(ih),'color',cs(ih,:))
        end
    else
        set(h,'color',cs(1,:))
    end
end
if not(ho)
    hold off
end
if nargout == 0
    clear h
end
function ytick(ticks,varargin)
set(gca,'ytick',ticks,varargin{:});

function yticklabel(labels)

set(gca,'yticklabel',labels)

function [varargout] = rep2struct(varargin)

% [s.target] = rep2struct(dat)
% replicate the value dat into each element of structure s in field target.
% if dat has same nb or elements as s, each element of dat goes into one
% element of s. if dat is more dimensional and doesn't have the same number
% of elements as s, and has same size along dimension 1, then pass each
% slice into s.target.

varargout = cell(1,nargout);
if numel(varargin) == 1
    dat = varargin{1};
    if numel(dat) == nargout
        for i = 1:nargout
            varargout{i} = dat(i);
        end
    elseif size(dat,1) == nargout
        for i = 1:nargout
            varargout{i} = dat(i,:);
        end
    else
        for i = 1:nargout
            varargout{i} = dat;
        end
    end
elseif numel(varargin) == nargout
    for i = 1:nargout
        varargout{i} = varargin{i};
    end
else
    error('Wrong number of arguments');
end
return


function s = setdef(s,d,keepempty)
% s = setdef(s,d)
% s = setdef(s,d,keepempty)
% Merges the two structures s and d recursively.
% Adding the default field values from d into s when not present or empty.
% Keeping order of fields same as in d
% if keepempty is provided and true, then empty fields in s will be left
% empty. otherwise they are populated with default values. default is
% false.
if not(exist('keepempty','var'))
    keepempty = 0;
end

if isstruct(s) && not(isempty(s))
    if not(isstruct(d))
        fields = [];
    else
        fields = fieldnames(d);
    end
    for i_f = 1:numel(fields)
        if isfield(s,fields{i_f})
            s.(fields{i_f}) = setdef(s.(fields{i_f}),d.(fields{i_f}),keepempty);
        else
            [s.(fields{i_f})] = d.(fields{i_f});
        end
    end
    if not(isempty(fields))
        fieldsorig = setdiff(fieldnames(s),fields);
        s = orderfields(s,[fields; fieldsorig]);
    end
elseif not(isempty(s)) || keepempty
    s = s;
elseif isempty(s)
    s = d;    
end
