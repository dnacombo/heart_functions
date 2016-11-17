function [beats, beats_time] = heart_peak_detect(cfg,data)

% [beats, beats_time] = heart_peak_detection(ECG,fs)
%
% Detect R peaks in ECG data sampled at fs Hz.
% Inputs:
%       ECG         vector of ECG data
%       fs          sampling frequency
% Outputs:
%       beats       samples where beats have been detected.
%       beats_time  time points where beats have been detected (assuming
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
% "data", then that should contain data that was already read from file in
% a previous call to FT_PREPROCESSING. In that case only the configuration
% options below apply.
%
% 
% The options are specified with
%   - Preprocessing options:
%     cfg.downsample      = 'no' or 'yes'  downsample the data to a specified rate (default = 'yes') 
%     cfg.downrate        = down-sampling rate (default = 300)
%     cfg.hpfilter        = 'no' or 'yes'  highpass filter (default = 'yes')
%     cfg.hpfreq          = highpass frequency in Hz (default = 1);
%
%   - Algorithm options (see description below):
%     cfg.thresh          = z-threshold for 1st step detection of R-peaks (default = 10)
%     cfg.mindist         = minimum inter beat interval in seconds (IBI) (default = 0.35)
%     cfg.corthresh       = proportion of maximum correlation (default = 0.6)
%
%   - Plotting options:
%     cfg.plotthresh      = open a figure to modify the threshold (default = 'no')
%     cfg.plotbeat        = open a figure to show the average ECG around R peak (default = 'no')
%     cfg.plotcorr        = open a figure to show the correlation level along the recording and enable editing threshold correlation (default = 'no')
%     cfg.plotall         = whether to do all of the above (default = 'no')
%
%
%
%   Algorithm:
%
%   The signal is first resampled and high pass filtered. The square of the
%   z-scored ECG is computed and a first detection is performed as the
%   peaks passing the cfg.thresh. Not all R peaks need to be selected at
%   this step. Just enough to create a template heart beat ECG (HB) is
%   necessary.
%   If cfg.plotthresh is true, then a figure is shown, allowing the user to
%   edit the threshold.
%   Then a template HB is computed (shown in a figure if cfg.plotbeat) and
%   convolved with the whole ECG time series. The resulting convolution is
%   normalized to have a maximum of 1 and beats are taken as peaks above
%   cfg.corthresh.
%   In both steps, a minimum distance between beats of cfg.mindist is
%   enforced.
%   
%  

% v0 Maximilien Chaumon November 2016
% based on previously undocumented anonymous version

narginchk(1,2)
structoutput = 0;
if nargin == 2 && isnumeric(cfg) && isnumeric(data)
    if not(isvector(cfg))
        error('ECG should be a one channel vector of data');
    end
    fs = data;
    data = [];
    data.fsample = fs;
    data.label = {'ECG'};
    data.trial = {cfg(:)'};
    data.time = {linspace(0,size(data.trial{1},2)/cfg.fsample,size(data.trial{1},2))};
    cfg = [];
elseif nargin == 1 && isstruct(cfg) && isstruct(ft_checkopt(cfg,'dataset','string'))
    data = ft_preprocessing(cfg);
elseif nargin == 2 && isnumeric(cfg) && isstruct(data) && isstruct(ft_checkopt(data,'fsample','double'))
    tmp = cfg;cfg = data;
    if not(isvector(tmp))
        error('ECG should be a one channel vector of data');
    end
    data = [];
    data.fsample = cfg.fsample;
    data.label = {'ECG'};
    data.trial = {tmp(:)'};
    data.time = {linspace(0,size(data.trial{1},2)/cfg.fsample,size(data.trial{1},2))};
elseif nargin == 2 && isstruct(cfg) && isstruct(data) && isstruct(ft_checkdata(data,'raw'))
    structoutput = 1;
else
    error('Wrong input');
end
if nargout == 2
    structoutput = 0;
end

ECG = data.trial{1};
if not(isvector(ECG))
    error('ECG input should be a vector or channel should be specified')
end


def = [];
def.downsample      = 'yes';% downsample?
def.downrate        = 300;  % down-sampling rate
def.hpfilter        = 'yes';
def.hpfreq          = 1;    % low bound of high pass filter of ECG
def.thresh          = 10;    % z-threshold for 1st step detection of R-peaks
def.mindist         = 0.35; % minimum IBI
def.corthresh       = 0.6;  % proportion of maximum correlation
def.plotall         = 0;
def.plotthresh      = 0;
def.plotbeat        = 0;
def.plotcorr        = 0;

cfg = setdef(cfg,def);



fsample_ori = data.fsample;

if istrue(cfg.downsample)
    tmp = [];
    tmp.resamplefs = cfg.downrate;
    tmp.detrend = 'yes';
    tmp.demean = 'yes';
    [data] = ft_resampledata(tmp,data);
    cfg.fsample = data.fsample;
end


tmp             = [];
tmp.hpfilter    = cfg.hpfilter;
tmp.hpfreq      = cfg.hpfreq;
data            = ft_preprocessing(tmp,data);

%% Here we use a threshold  z-score method to find the R-peaks

% standardize ecg
ECG = data.trial{1};
ECG2z = nanzscore(ECG).^2;

thresh = cfg.thresh; % default z-threshold
[beats, v] = peakdetect2(ECG2z, thresh, cfg.fsample .* cfg.mindist);
while istrue(cfg.plotthresh) || istrue(cfg.plotall)
    figure(47894);clf
    plot(ECG2z);
    hold on;
    scatter(beats,v);
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
        [beats, v] = peakdetect2(ECG2z, thresh, cfg.fsample .* cfg.mindist);
    end
end



%% We now build the template and compute the correlation with the ecg channel

% build template based on identified r-peaks
HB_bound = round(.5 * cfg.fsample);
HB = NaN(numel(beats),2*HB_bound+1);
for ii=1:length(beats)
    if beats(ii)-HB_bound > 1 && beats(ii)+HB_bound < length(ECG)
        HB(ii,:) = ECG(beats(ii)-HB_bound:beats(ii)+HB_bound);
    end
end

mHB = mean(HB,1);

% we'll assume that signal at dt seconds before detected R-peaks should be
% lower than R-peaks. 50 ms seems reasonable.
dt = round(.05 * cfg.fsample);
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
[beats, v] = peakdetect2(cr(1001:end-1000), thresh, cfg.fsample .* cfg.mindist);

while istrue(cfg.plotcorr) || istrue(cfg.plotall)
    IBI_s = diff(beats)/cfg.fsample;
    plotIBI(ecg_n,cr,beats,v,thresh,IBI_s);
    s1 = questdlg('Would you like to change the threshold?');
    switch s1
        case 'Yes'
            [dum,thresh] = ginput(1);
            [beats, v] = peakdetect2(cr(1001:end-1000), thresh, cfg.fsample .* cfg.mindist);
        case 'No'
            break
        case 'Cancel'
            error('Don''t click cancel unless you want to cancel...')
    end
end
while istrue(cfg.plotcorr) || istrue(cfg.plotall)
    plotIBI(ecg_n,cr,beats,v,thresh,IBI_s);
    
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
            for ii = beats(IBI_out)
                while istrue(cfg.plotcorr) || istrue(cfg.plotall)
                    h = figure(478490);clf;
                    plot(ecg_n);
                    hold on
                    %                     plot(cr(1001:end-1000),'r')
                    scatter(beats,ecg_n(beats))
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
                        beats(end+1) = closest; beats = sort(beats);
                    else
                        % find closest detected peak
                        [dum,closest] = min(abs(beats-x));
                        beats(closest) = [];
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
close
if istrue(cfg.downsample)
    beats = round(beats*fsample_ori/cfg.fsample);
end

beats_time = beats / fsample_ori;

if structoutput
    beats = struct('R_sample',beats,'R_time',beats_time);
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



function [p, v] = peakdetect2(dat, val, mindist)

% PEAKDETECT2 detects peaks above a certain threshold in single-channel data
%
% Use as
%   [pindx, pval] = peakdetect(signal, min, mindist)
%
% mindist is optional, default is 1
%
% See also PEAKDETECT, PEAKDETECT3

% Copyright (C) 2000, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: peakdetect2.m 952 2010-04-21 18:29:51Z roboos $



% if nargin<3
%   mindist=1;
% end

i = find(dat>val);
m = dat(i);
d = diff(i);
jump = (d>mindist);
p = [];

sect=1;
while sect<=length(d)
  if jump(sect)
    p = [p i(sect)];
  else
    s = [];
    while ~jump(sect) & sect<length(d)
      s = [s sect];
      sect = sect + 1;
    end
    [lm, li] = max(m(s));
    p = [p i(s(li))];
  end
  sect = sect+1;
end

if nargout>1
  v = dat(p);
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



