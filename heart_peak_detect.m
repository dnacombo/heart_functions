function [HeartBeats, R_time] = heart_peak_detect(cfg,data)

% [R_sample, R_time] = heart_peak_detection(ECG,fs)
%
% Detect R peaks in ECG data sampled at fs Hz.
% Inputs:
%       ECG         vector of ECG data
%       fs          sampling frequency
% Outputs:
%       R_sample    samples where beats have been detected.
%       R_time      time points where beats have been detected (assuming
%                   time starts at 0)
%
% [HeartBeats] = heart_peak_detect(cfg)
% [HeartBeats] = heart_peak_detect(cfg,data)
%
% If you are calling heart_peak_detect with only the configuration as first
% input argument and the data still has to be read from file, you should
% specify
%   cfg.dataset      = string with the filename
%   cfg.channel      = the one channel to be read and taken as ECG channel
%
% If you are calling heart_peak_detect with also the second input argument
% "data", then that should either
%   - be a vector of ECG data points. first input cfg.fsample should
%   contain sampling frequency in this case.
%   - contain data that was read from file in a previous call to
%   FT_PREPROCESSING.
%
% In either case, the configuration options below apply.
%
% 
% The options are specified with
%   - Preprocessing options:
%     cfg.hpfilter        = 'no' or 'yes'  highpass filter (default = 'yes')
%     cfg.hpfilttpye      = filter type (see ft_preprocessing, default = 'firws')
%     cfg.hpfreq          = highpass frequency in Hz (default = 1);
%     cfg.lpfilter        = 'no' or 'yes'  lowpass filter (default = 'yes')
%     cfg.lpfilttpye      = filter type (see ft_preprocessing, default = 'firws')
%     cfg.lpfreq          = lowpass frequency in Hz (default = 100);
%
%   - Algorithm options (see description below):
%     cfg.thresh          = z-threshold for 1st step detection of R-peaks (default = 10)
%     cfg.mindist         = minimum inter beat interval in seconds (IBI) (default = 0.35)
%     cfg.corthresh       = proportion of maximum correlation (default = 0.6)
%     cfg.QRmax           = maximum duration between Q and R peaks in seconds (default = 0.05)
%     cfg.RSmax           = maxumum duration between R and S peaks in seconds (default = 0.1)
%     cfg.QTmax           = maximum QT interval in seconds (default = 0.42)
%
%   - Plotting options:
%     cfg.plotthresh      = open a figure to modify the threshold (default = 'no')
%     cfg.plotbeat        = open a figure to show the average ECG around R peak (default = 'no')
%     cfg.plotcorr        = open a figure to show the correlation level along the recording and enable editing threshold correlation (default = 'no')
%     cfg.plotfinal       = open a figure to show final results with all peaks found (default = 'no')
%     cfg.plotall         = whether to do all of the above (default = 'no')
%
%   Algorithm:
%
%   The signal is first high and low pass filtered (default 1-100 Hz). The
%   square of the z-scored ECG is computed and a first detection is
%   performed as the peaks passing the cfg.thresh. Not all R peaks need to
%   be selected at this step. Just enough to create a template heart beat
%   ECG (HB) is necessary.
%   If cfg.plotthresh is true, then a figure is shown, allowing the user to
%   edit the threshold.
%   Then a template HB is computed (shown in a figure if cfg.plotbeat) and
%   convolved with the whole ECG time series. The resulting convolution is
%   normalized to have a maximum of 1 and beats are taken as peaks above
%   cfg.corthresh.
%   In both steps, a minimum distance between beats of cfg.mindist is
%   enforced.
%   Other peaks are found based on each R peak. Q is the minimum within 50
%   ms before R, S is the minimum within 100 ms after R, and T is the
%   maximum between the S peak and a maximum QT interval of 420 ms (a rough
%   standard...).
%  

% v0 Maximilien Chaumon November 2016
% based on previously undocumented anonymous version

narginchk(1,2)
if isnumeric(cfg) % first input method
    if not(isvector(cfg))
        error('ECG should be a one channel vector of data');
    end
    if not(exist('data','var'))
        error('Please provide sampling rate');
    end
    fs = data;
    data = [];
    data.fsample = fs;
    data.label = {'ECG'};
    data.trial = {cfg(:)'};
    data.time = {linspace(0,size(data.trial{1},2)/fs,size(data.trial{1},2))};
    cfg = [];
    cfg.structouput = 0;
elseif isstruct(cfg) % second input method
    if nargin == 1 
        data = ft_preprocessing(cfg);
    elseif nargin == 2
        if isnumeric(data)% data input as vector
            if ~isvector(data)
                error('ECG should be a one channel vector of data');
            end
            if ~isstruct(ft_checkopt(cfg,'fsample','double'))
                error('Please provide sampling frequency in cfg.fsample')
            end
            tmp = data;
            data = [];
            data.fsample = cfg.fsample;
            data.label = {'ECG'};
            data.trial = {tmp(:)'};
            data.time = {linspace(0,size(data.trial{1},2)/data.fsample,size(data.trial{1},2))};
            cfg = rmfield(cfg,'fsample');
        elseif isstruct(ft_checkdata(data,'datatype','raw'))
            if isfield(cfg,'channel')
                data = ft_preprocessing(cfg,data);
            end
            if not(isvector(data.trial{1}))
                error('Channel should be specified')
            end
        end
    end
else
    error('Bad input');
end

data_in = data;

def = [];
def.hpfilter        = 'yes';
def.hpfreq          = 1;    % low bound of high pass filter of ECG
def.hpfilttype      = 'firws';
def.lpfilter        = 'yes';
def.lpfilttype      = 'firws';
def.lpfreq          = 100;
def.thresh          = 10;    % z-threshold for 1st step detection of R-peaks
def.mindist         = 0.35; % minimum IBI
def.corthresh       = 0.6;  % proportion of maximum correlation
def.QRmax           = 0.05;  
def.RSmax           = 0.1;
def.QTmax           = 0.42;
def.structoutput    = 1; 
def.plotall         = 0;
def.plotthresh      = 0;
def.plotbeat        = 0;
def.plotcorr        = 0;
def.plotfinal       = 0;

cfg = setdef(cfg,def);

if cfg.lpfreq > data.fsample/2
    error(['Lowpass filter frequency too high. Set cfg.lpfreq below ' num2str(data.fsample/2)])
end
tmp             = [];
tmp.hpfilter    = cfg.hpfilter;
tmp.hpfreq      = cfg.hpfreq;
tmp.hpfilttype  = cfg.hpfilttype;
tmp.lpfilter    = cfg.lpfilter;
tmp.lpfreq      = cfg.lpfreq;
tmp.lpfilttype  = cfg.lpfilttype;
data            = ft_preprocessing(tmp,data);

ECG = data.trial{1};
time = data.time{1};

%% find R peaks

% standardize ecg
ECG2z = nanzscore(ECG).^2;

thresh = cfg.thresh; % default z-threshold
[R_sample, R_value] = peakseek(ECG2z, thresh, data.fsample .* cfg.mindist);
while istrue(cfg.plotthresh) || istrue(cfg.plotall)
    figure(47894);clf
    plot(ECG2z);
    hold on;
    scatter(R_sample,R_value);
    hline(thresh,':r');
    xlabel('samples');
    ylabel('zscore');
    zoom;
    
    title(sprintf(['Creating a template HB\nWe should have enough HB (don''t need all)\nclick to change threshold\nright-click to confirm'], thresh));
    
    [x,y,but] = ginput(1);
    if ~isempty(but) && but > 1
        break
    else
        thresh = y;
        [R_sample, R_value] = peakseek(ECG2z, thresh, data.fsample .* cfg.mindist);
    end
end


%% We now build the template and compute the correlation with the ecg channel

% build template based on identified r-peaks
HB_bound = round(.5 * data.fsample);
HB = NaN(numel(R_sample),2*HB_bound+1);
for ii=1:length(R_sample)
    if R_sample(ii)-HB_bound > 1 && R_sample(ii)+HB_bound < length(ECG)
        HB(ii,:) = ECG(R_sample(ii)-HB_bound:R_sample(ii)+HB_bound);
    end
end

mHB = nanmean(HB,1);

% we'll assume that signal at dt seconds before detected R-peaks should be
% lower than R-peaks. dt = 50 ms seems reasonable.
dt = round(.05 * data.fsample);
if ~(mean(mHB(HB_bound-5:HB_bound+5)) > mean(mHB(HB_bound - dt:HB_bound - dt + 10)))
    % flip
    ECG = -ECG;
    mHB = -mHB;
    HB = -HB;
end
while istrue(cfg.plotbeat) || istrue(cfg.plotall)
    figure(47894);clf
    t = linspace(-HB_bound,HB_bound,numel(mHB));
    hall = plot(t,HB','color',[.8 .8 .8]);
    hold on
    hm = plot(t,mHB);
    vline(0,'r');
    legend([hall(1),hm],{'individual ECG','mean ECG'},'location','best')
    axis tight
    title('r peak points up ? ''left mouse'' = no ''right mouse'' = yes');
    [x,y,but] = ginput(1);
    if ~isempty(but) && but > 1
        break
    else
        ECG = -ECG;
        mHB = -mHB;
        HB = -HB;
    end
end

ecg_n = ECG./max(ECG); % normalized ecg

% pad ecg
ecg_pad = [zeros(1,1000) ECG zeros(1,1000)];
cr = zeros(size(ecg_pad));

% compute correlation
for ii=1:length(ecg_pad)-length(mHB)
    cr(ii+round(length(mHB)/2)-1) = sum(ecg_pad(ii:ii+length(mHB)-1).*mHB);
end
cr = cr./max(cr); % normalize correlation to 1

% find peaks in correlation
thresh = cfg.corthresh;
[R_sample, R_value] = peakseek(cr(1001:end-1000), thresh, data.fsample .* cfg.mindist);

while istrue(cfg.plotcorr) || istrue(cfg.plotall)
    IBI_s = diff(R_sample)/data.fsample;
    plotIBI(ecg_n,cr,R_sample,R_value,thresh,IBI_s);
    s1 = questdlg('Would you like to change the threshold?');
    switch s1
        case 'Yes'
            [dum,thresh] = ginput(1);
            [R_sample, R_value] = peakseek(cr(1001:end-1000), thresh, data.fsample .* cfg.mindist);
        case 'No'
            break
        case 'Cancel'
            error('Don''t click cancel unless you want to cancel...')
    end
end
while istrue(cfg.plotcorr) || istrue(cfg.plotall)
    plotIBI(ecg_n,cr,R_sample,R_value,thresh,IBI_s);
    
    s2 = questdlg('Are there still outliers present? ');
    switch s2
        case 'Yes'
            axes(findobj(gcf,'tag','hhist'))
            title('click minimal bound of correct IBI')
            [lw,~] = ginput(1);
            title('click maximal bound of correct IBI')
            [hg,~] = ginput(1);
            % find outlier peaks
            IBI_out = IBI_s < lw | IBI_s > hg;
            for ii = R_sample(IBI_out)
                while istrue(cfg.plotcorr) || istrue(cfg.plotall)
                    h = figure(478490);clf;
                    plot(ecg_n);
                    hold on
                    %                     plot(cr(1001:end-1000),'r')
                    scatter(R_sample,ecg_n(R_sample))
                    xlim([ii-500 ii+500]);
                    xlabel('samples')
                    ylabel('a.u.')
                    title('Add with left mouse. Remove with right mouse. Move on with Escape.')
                    [x,y,but] = ginput(1);
                    x = round(x);
                    if but == 27
                        break
                    elseif but == 1
                        % find closest maximum
                        X = [1:numel(ecg_n);ecg_n*numel(ecg_n)]';
                        [closest] = dsearchn(X,[x y*numel(ecg_n)]);
                        scatter(closest,ecg_n(closest))
                        R_sample(end+1) = closest; R_sample = sort(R_sample);
                    else
                        % find closest detected peak
                        [dum,closest] = min(abs(R_sample-x));
                        R_sample(closest) = [];
                    end
                end
            end
            try
                close(h)
            end
        otherwise
            break
    end
end

if cfg.structoutput
    %% find Q S T
    Q_sample = NaN(size(R_sample));S_sample = NaN(size(R_sample));T_sample = NaN(size(R_sample));
    for i_R = 1:numel(R_sample)
        bnds = max(1,round(R_sample(i_R) - cfg.QRmax * data.fsample)):R_sample(i_R);
        ECGtmp = ECG(bnds);
        [v,p] = min(ECGtmp);
        Q_sample(i_R) = p + bnds(1) - 1;
        bnds = R_sample(i_R):min(numel(ECG),round(R_sample(i_R) +  cfg.RSmax * data.fsample));
        ECGtmp = ECG(bnds);
        [v,p] = min(ECGtmp);
        S_sample(i_R) = p + bnds(1) - 1;
        QTmax = 0.42;
        bnds = S_sample(i_R):min(numel(ECG),round(Q_sample(i_R)+cfg.QTmax*data.fsample));
        ECGtmp = ECG(bnds);
        [v,p] = max(ECGtmp);
        T_sample(i_R) = p + bnds(1) - 1;
    end
    Q_time = skipnans(data_in.time{1},Q_sample);
    R_time = skipnans(data_in.time{1},R_sample);
    S_time = skipnans(data_in.time{1},S_sample);
    T_time = skipnans(data_in.time{1},T_sample);
    
    [HeartBeats.Q_sample] = rep2struct(Q_sample);
    [HeartBeats.Q_time] = rep2struct(Q_time);
    [HeartBeats.R_sample] = rep2struct(R_sample);
    [HeartBeats.R_time] = rep2struct(R_time);
    [HeartBeats.S_sample] = rep2struct(S_sample);
    [HeartBeats.S_time] = rep2struct(S_time);
    [HeartBeats.T_sample] = rep2struct(T_sample);
    [HeartBeats.T_time] = rep2struct(T_time);
    if istrue(cfg.plotfinal) || istrue(cfg.plotall)
        figure;
        hold on
        set(gca,'colororder',[0         0    1.0000
            0    0.5000         0
            0.8000         0         0
            0    0.7500    0.7500
            0.7500         0    0.7500]);
        
        plot(data_in.time{1},data_in.trial{1},'k');
        
        todel = isnan(Q_sample);
        toplot = [Q_time(~todel); data_in.trial{1}(Q_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(R_sample);
        toplot = [R_time(~todel); data_in.trial{1}(R_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(S_sample);
        toplot = [S_time(~todel); data_in.trial{1}(S_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        todel = isnan(T_sample);
        toplot = [T_time(~todel); data_in.trial{1}(T_sample(~todel))];
        plot(toplot(1,:),toplot(2,:),'.','markersize',10);
        legend({'ECG','Q','R','S','T'});
    end
    clear beats_time
else
    HeartBeats = R_sample;
end

function [b] = skipnans(a,idx)

b = NaN(size(idx));
for i = 1:numel(idx)
    if ~isnan(idx(i))
        b(i) = a(idx(i));
    end
end
        

function plotIBI(ecg_n,cr,p,v,thresh,IBI_s)
figure(47894);clf
wholescreen = get(0,'ScreenSize');
pp = wholescreen;
pp(2) = wholescreen(4)/2;
pp(4) = wholescreen(4)/2;
set(gcf,'position',pp)

subplot(2,1,1)
plot(ecg_n,'b');
hold on
plot(cr(1001:end-1000),'r')
scatter(p,v)

axis tight
hline(thresh,'linestyle','--','color','k')
xlabel('samples')
ylabel('a.u.')
legend('ECG','corr','R-peak','location','eastoutside')
zoom
title('zoom in to check signal')

hhist = subplot(2,2,3);
hist(IBI_s,30)
xlabel('IBI (s)')
set(hhist,'tag','hhist')


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


