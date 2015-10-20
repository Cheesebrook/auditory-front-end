%plotbinCueHistc.m
%
% plotbinCueHistc(Nhist,gravityFunction,varargin)
%
%This function plots a histogram 'Nhist' calculated by function
%binCueHistc.m 
%Either a surf or image plot can be plotted specified by 'method'.
%
%inputs:
% Nhist - The histogram with dimensions [nbins x Nchan] calculated by binCueHistc.m 
%        The plot will be along dimensions [x y Ndist]
% optional:
% 'method' - 'IAMGE' uses imagesc, 'SURF' uses surf plot
% 'gravity' - bool: plot gravity of histogram(1)
% 'clim' - color scale limits for IMAGE plot (see imagesc.m)
% 'gravlim' - ylim for gravity function plot
% 'brel' - bool: relative occurance (normalise data) (1)
% 'bgravity' - bool: relative occurance (normalise data) (1)
% 'ymax' - maximum value for y-vector
% 'labels' - xlabel and ylabel
% 'S' - plot settings: FontSize
%
%-----------------------------------------------
%Revision:
% 17-10-2013 by JK: including flexible inputs via varargin 
% 24-03-2014 by JK: - changed 'thresh' to 'actbands' (in order to omit a second estimation of the activated bands inside this function)
%                   - changed histogram-calculation from 'hist' to 'histc' 
%-----------------------------------------------
%
%by Johannes K?sbach (JK), DTU, CAHR, 08-October-2013
%
% see also: binCueHistc.m, surf.m, imagesc.m

function plotbinCueHistc(Nhist,gravityFunction,varargin)

%% Definitions (fixed Parameters)
gravplwidth = 0.1; %width of plot for gravity
shift = 0.3; %shift in width for imagesc and pos(1) for gravity plot
xtick = 5; %every 'xtick' tick on xaxis
% fsax = 14; %FontSize Axis -> substituted by S.fsax

%% check inputs
if nargin<2||isempty(gravityFunction) 
    % consider only activated bands
    [~, Nchanact] = getActBands(Nhist);
    
    %generate gravity Function
    gravityFunction = nansum(Nhist,2)/Nchanact;
end

r = struct(varargin{:});
try 
    clim = r.clim;
catch
    % clim = [0 1]; %limits for colorbar in imagesc
    clim = [0 0.3]; %limits for colorbar in imagesc for FA2014 presentation
end
try
    gravlim = r.gravlim;
catch
    gravlim = 0.2; %limit of x-axis (y-axis before rotation on view) of plot for gravity
end
try
    color_flag = r.color_flag;
catch
    color_flag = 0; %bool: plot colorbar(1) 
end
try
    method  = r.method;
catch
    method = 'IMAGE'; %imagesc plot
end
try
    brel  = r.brel;
catch
    brel = 1; %relative occurance (normalise data)
end
try
    bgravity  = r.bgravity;
catch
    bgravity = 1; %relative occurance (normalise data)
end
try
    ymax = r.ymax;
catch
    ymax = max(max(abs(Nhist))); %ylim according to data
end
try
    labels = r.labels;
catch
    labels = {'Filterbank Channel';'Binaural Cue'}; %labels
end
try
    xticklabel = r.xticklabel; %use for cfHz
catch
    %do nothing
end
try
    markers = r.markers; %vertical lines to mark a certain frequncy region
catch
    markers = 0;
end
try
    PLOT = r.PLOT;
catch
    load S %common plot settings
    PLOT = S;
end
% end

if length(labels)<3
    if(brel)
        labels{3} = 'Rel. occurance'; %(normalised data)
    else
        labels{3} = 'Abs. occurance';
    end
end

%% Plots

%number of containers 'nbins' and frequency channels 'Nchan'
nbins = size(Nhist,1);
Nchan = size(Nhist,2);

%x and y coordinate
x = 1:Nchan;
y = linspace(-ymax,ymax,nbins);

figure
if bgravity
    subplot(1,2,1)
end
switch method
    case 'SURF'
        [Xmesh, Ymesh] = meshgrid(x,y); %meshgrid
        surf(Xmesh,Ymesh,Nhist);
        zlabel(labels{3},'FontSize',PLOT.fsztxt)
        plotone = gca;
    case 'IMAGE'
        imagesc(x,y,Nhist,clim);
        colormap('bone')
        colormap(flipud(colormap))
        if color_flag
            colorbar
        end
        if markers
           hold on
           plot(markers(1),y,'Color','k','Linewidth',PLOT.lw)
           plot(markers(2),y,'Color','k','Linewidth',PLOT.lw)
           hold off
        end
        title(labels{3},'FontSize',PLOT.fsztxt)
        plotone = gca;
end

%labels
xlabel(labels{1},'FontSize',PLOT.fsztxt)
ylabel(labels{2},'FontSize',PLOT.fsztxt)

%axes
xlim([1 Nchan])
ylim([-ymax ymax])
if brel
    zlim(clim)
end

% if ~exist('xticklabel','var')
%     xticklabel = x;
% end
% set(plotone,'XTick',xtick:xtick:length(xticklabel))
% set(plotone,'XTickLabel',xticklabel(xtick:xtick:end))

% Managing frequency axis ticks for auditory filterbank
if exist('xticklabel','var')
    cfHz = xticklabel;
    % Find position of y-axis ticks
    M = size(cfHz,2);  % Number of channels
    n_points = 500;    % Number of points in the interpolation
    interpolate_ticks = spline(1:M,cfHz,linspace(0.5,M+0.5,n_points));

    % Restrain ticks to signal range (+/- a half channel)
    aud_ticks = [100 250 500 1000 2000 4000 8000 16000 32000];
    aud_ticks = aud_ticks(aud_ticks<=interpolate_ticks(end));
    aud_ticks = aud_ticks(aud_ticks>=interpolate_ticks(1));
    n_ticks = size(aud_ticks,2);        % Number of ticks
    ticks_pos = zeros(size(aud_ticks)); % Tick position

    % Find index for each tick
    for ii = 1:n_ticks
        jj = find(interpolate_ticks>=aud_ticks(ii),1);
        ticks_pos(ii) = jj*M/n_points;
    end
    
    % Set up y-axis
    set(plotone,'XTick',ticks_pos,...
        'XTickLabel',aud_ticks,'fontsize',PLOT.fszax) %,...
        %'fontname',p.map('ftype'))
        
else
    xticklabel = x;
    set(plotone,'XTick',xtick:xtick:length(xticklabel))
    set(plotone,'XTickLabel',xticklabel(xtick:xtick:end))
    set(gca,'FontSize',PLOT.fszax)
end

%gravity
if bgravity
    subplot(1,2,2)
    plot(y,gravityFunction,'LineWidth',PLOT.lw,'Color','k')
    view(90,90)
%     title('Centre of gravity','FontSize',S.fstxt)
    plottwo = gca;
    pos1 = get(plotone,'Position');
    pos2 = get(plottwo,'Position');
    hpytick = get(plotone,'YTick');
    set(plotone,'Position',[pos1(1) pos1(2) pos1(3)+shift pos1(4)])
    set(plottwo,'Position',[pos2(1)+shift-0.08 pos2(2) gravplwidth pos2(4)])
    set(gca,'XTick',hpytick)
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',{'0';'';num2str(gravlim)})
    xlim([-ymax ymax])
    ylim([0 gravlim])
    set(gca,'FontSize',PLOT.fszax)
end
        
end
%EOF