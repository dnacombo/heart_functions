function HeartBeats = heart_SV_calc(cfg,data,HeartBeats)

% [HeartBeats] = heart_SV_calc(cfg)
% [HeartBeats] = heart_SV_calc(cfg,data)
% [HeartBeats] = heart_SV_calc(cfg,data,HeartBeats)
%
%   The results are returned in HeartBeats, a vector structure with as many
%   elements as detected R peaks in the ECG. Retrieve Stroke volume for
%   each heart beat with : [HeartBeats.SV]
% 
% If you are calling heart_SV_calc with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.dataset      = string with the filename
%
% If you are calling heart_SV_calc with also the second input argument
% "data", then that should contain data that was read from file in a
% previous call to FT_PREPROCESSING. The configuration options below apply.
%
% If HeartBeats is provided as a third argument, it is assumed to be the
% output of heart_peak_detect. No detection is computed on the ECG.
%
% The options are specified with:
%     - ECG options:
%         cfg.ECG.channel     = channel to use for ECG (default = 'ECG1')
%         cfg.ECG.scale       = scale, to display along with Z0 (default = -4000)
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
%         cfg.LVETmax         = maximal left ventricular ejection time in seconds (default = .4)
%         cfg.LVE_startfracmax_dZdt = fraction of maximum dZdt value to use to detect ejection start (default = 0.15, Kubicek et al. 1970)
%         cfg.rho             = resistivity of blood in Ohms (default = 135)
%         cfg.L               = distance between the two middle electrodes in cm (default = 30)
% 
%     - Ploting options:
%         cfg.plotallmeas     = whether or not to plot all measures (default = 'no')
%         cfg.interactive     = whether or not an interactive plot should be opened at the end to allow changint outlier peaks (default = 'yes')
%
%
%   Algorithm:
%
% Peaks are first detected on the ECG signal using heart_peak_detect with
% the options specified. The last heart beat is discarded to ensure all
% heart cycles are complete. The impedance channel is low pass filtered,
% dZ/dt is computed.
% For each detected heart cycle, start with the Q point. Find the maximum
% dZ/dt (dZdt_max) following this Q point (within QCmax seconds). Find the time at
% which dZ/dt reaches LVE_startfracmax_dZdt * dZdt_max. This is the B
% point (start of the left ventricular ejection (LVE)). Find the time
% following dZdt_max where dZ/dt reaches a minimum value (within LVETmax
% seconds of the B point). This is the X point (end of LVE). Compute LVET,
% the duration in seconds of LVE.
% Compute Z, the mean Z0 during LVE.
% Compute Stroke volume with the following formula:
%
% SV = rho .* L .^ 2 ./ Z .^ 2 .* LVET .* dZdt_max
%

% v0 Maximilien Chaumon November 2016


def             = [];
def.ECG.channel = 'ECG1';
def.ECG.scale   = -4000;
def.ECG.show    = 'yes';

def.ECG_peakcfg    = [];

def.Z0.channel  = 'EBI_-_Magnitude';
def.Z0.lpfilter = 'yes';
def.Z0.lpfreq   = 15;
def.Z0.lpfilttype = 'but';
def.Z0.factor   = 5;
def.Z0.show     = 'yes';

def.sound.channel = 'HeartSound';
def.sound.bpfilter = 'yes';
def.sound.bpfreq = [40 60];
def.sound.tubelength = 152;

def.QCmax       = .3;
def.LVETmax     = .4;
def.LVE_startfracmax_dZdt = 0.15;
def.rho         = 135;
def.L           = 30;
def.SV          = @(rho,L,Z,LVET,dZdt_max_after_HB) rho .* L .^ 2 ./ Z .^ 2 .* LVET .* dZdt_max_after_HB;

def.plotallmeas = 'no';
def.interactive = 'yes';

cfg = setdef(cfg,def,1);

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


if hassound
    %% detect heart sounds (S1, S2)
    dataSound = ft_preprocessing(cfg.sound,data);
    HeartSound_bp = dataSound.trial{1};
    soundlag =  round(cfg.sound.tubelength/33000 * data.fsample);% in samples
    HeartSound_bp = [HeartSound_bp(soundlag:end) NaN(1,soundlag-1)];

    HeartSoundz = nanzscore(HeartSound_bp).^2;
    HeartSoundz(lastHB.R_sample:end) = NaN;
    % SoundPks = peakdetect2(HeartSoundz,1,.1*data.fsample);
    % NB:   peakdetect2 detects peak peak (highest peak above thresh)
    %       peakdetect3 detects onset (first peak above thresh)
    % may want to refine this...
    SoundPks_idx = peakseek(HeartSoundz,1,.1*data.fsample);
    SoundPks_time = time(SoundPks_idx);

end

%% Impedance: filtering
dataZ0 = ft_preprocessing(cfg.Z0,data);
Z_0 = - dataZ0.trial{1} * cfg.Z0.factor; % use - by convention
Z_0(lastHB.R_sample:end) = NaN;


%% Impedance: dZ/dt
dZdt = diff(Z_0);
dZdt = [dZdt NaN]; % pad with one NaN
dZdt = dZdt * data.fsample; % /dt
dZdt(lastHB.R_sample:end) = NaN;

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
    LVET_time(i_R) = time(LVET_sample(i_R));
    Z(i_R) = mean(Z_0(B_sample(i_R):X_sample(i_R)),2);
    
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

[HeartBeats.SV] = rep2struct(cfg.SV(cfg.rho,cfg.L,[HeartBeats.Z],[HeartBeats.LVET_time],[HeartBeats.dZdt_max_value]));
HeartBeats(1).cfg = cfg;

if istrue(cfg.plotallmeas)
    hheart = plot(time,cfg.ECG.scale*ECG);
    vline([HeartBeats.R_time],':','color',get(hheart,'color'))
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
    GD.ECG = cfg.ECG.scale*ECG;
    GD.Z_0 = Z_0;
    GD.dZdt = dZdt;
    GD.HeartSoundz = HeartSoundz;
    GD.HeartBeats = HeartBeats;
    GD.iout = 1;
    GD.LVE_ori = [GD.HeartBeats.B_time ;GD.HeartBeats.X_time]' - [GD.HeartBeats.R_time ;GD.HeartBeats.dZdt_max_time]';

    [ioutliers] = find_outliers(GD.LVE_ori);
    [GD.allout,GD.allout_type] = find(ioutliers);

    guidata(hfig,GD);

    plotoutliers(hfig,[]);
    uiwait(hfig);
end

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
        timearound = winaround + GD.HeartBeats(GD.allout(GD.iout)).B_time;
    case 2
        timearound = winaround + GD.HeartBeats(GD.allout(GD.iout)).X_time;
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

vline([HeartBeats.B_time],'color','k');
plot([HeartBeats.B_time],GD.dZdt([HeartBeats.B_sample]),'.k')

vline([HeartBeats.X_time],'color','k');
plot([HeartBeats.X_time],GD.dZdt([HeartBeats.X_sample]),'.k')

switch GD.allout_type(GD.iout)
    case 1
        vline(GD.HeartBeats(GD.allout(GD.iout)).B_time,'r')
    case 2
        vline(GD.HeartBeats(GD.allout(GD.iout)).X_time,'r')
end
xlim(timearound)
xlabel('Time (s)')

set(gca,'ButtonDownFcn',@editpeaks);

subplot(1,5,5)
LVE = [GD.HeartBeats.B_time ;GD.HeartBeats.X_time]' - [GD.HeartBeats.R_time ;GD.HeartBeats.dZdt_max_time]';
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
allpoints = [[GD.HeartBeats.B_time]' [GD.dZdt([GD.HeartBeats.B_sample])]';
    [GD.HeartBeats.X_time]' [GD.dZdt([GD.HeartBeats.X_sample])]'];
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
        GD.HeartBeats(iBeat).B_time = nudat(1);
        GD.HeartBeats(iBeat).B_sample = nupoint;
    case 2
        % LVE end
        GD.HeartBeats(iBeat).X_time = nudat(1);
        GD.HeartBeats(iBeat).X_sample = nupoint;
end

guidata(hObject,GD);
plotoutliers(hObject);


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
