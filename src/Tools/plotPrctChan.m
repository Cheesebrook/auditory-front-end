%plotPrctChan.m
%
% ------------------------ Modelling ASW -------------------------------- %
%
% This framework concerns modelling apparent source width (ASW) perception
% based on data measured and presented in Käsbach et al. (2013-2015)
%
% by Johannes Käsbach (JK), Centre for Applied Hearing Research (CAHR),
% Technical University of Denmark (DTU), johk@elektro.dtu.dk
%
% ----------------------------------------------------------------------- %
%
% plotPrctChan(prctChan,varargin)
%
% This function plots the percentiles per frequency channel 'prctChan' of 
% a binaural cue (itd, ild or ic) or another ASW-estimation calculated in
% calcASW.m
%
% inputs:
% prctChan - Percentiles per Channel [prct x Nchan]
%
% optional:
%
%
% --------History-----------
% - 30. September 2015: Modelling ASW framework
% - 16. Oktober 2015: 
% --------------------------
%
% see also: calcASW.m

function plotPrctChan(prctChan,varargin)

%% check inputs
r = struct(varargin{:});
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
    legendtag = r.legendtag;
catch
    legendtag = {'10 %','30 %','50 %','70 %','90 %'}; %lengend tags
end
try
    blegend = r.blegend;
catch
    blegend = 1; %bool: plot lengend
end
try
    xticklabel = r.xticklabel; %use for cfHz
catch
    %do nothing
end
try
    cfSel = r.cfSel; %selecting frequency bands for averaging of cohAv
catch
    cfSel = 1:size(prctChan,2);
end
try
    bvertmarkers = r.bvertmarkers; %vertical lines to mark a certain frequncy region
catch
    bvertmarkers = 0;
end
try
    choosePrct = r.choosePrct; %choose the percentiles to plot
    if choosePrct > size(prctChan,1)
        error('You cannot plot more percentiles than the array prctChan provides!')
    elseif choosePrct == 0;
        choosePrct = 1:size(prctChan,1);
    end
catch
    choosePrct = 1:size(prctChan,1);
end
try
    save_plot = r.save_plot; %bool: save plots(1)
catch
    save_plot = 0;
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

%select frequency channels
prctChan = prctChan(:,cfSel);

%number of containers 'nbins' and frequency channels 'Nchan'
Nprct = size(prctChan,1);
Nchan = size(prctChan,2);

%x vector
x = 1:Nchan;

figure
hold on
for ii = choosePrct
    plot(x,prctChan(ii,:),'Marker',PLOT.marker2{ii},'MarkerSize',PLOT.ms,...
        'Color',PLOT.colors(ii,:),'LineStyle','none')
end
plotone = gca;
if bvertmarkers
    hold on
    plot(repmat(x,2,1),repmat([-ymax ymax],length(x),1)','Color','k',...
        'Linewidth',1,'LineStyle',':')
    hold off
end

%labels
if blegend
    legend(legendtag,'Location','SouthEast')
end
xlabel(labels{1},'FontSize',PLOT.fsztxt)
ylabel(labels{2},'FontSize',PLOT.fsztxt)
title(labels{3},'FontSize',PLOT.fsztxt)

%axes
xlim([1 Nchan])
ylim([-ymax ymax])

% Managing frequency axis ticks for auditory filterbank
if exist('xticklabel','var')
    cfHz = xticklabel(:,cfSel)/1e3;
    % Find position of y-axis ticks
    M = size(cfHz,2);  % Number of channels
    n_points = 500;    % Number of points in the interpolation
    interpolate_ticks = spline(1:M,cfHz,linspace(0.5,M+0.5,n_points));

    % Restrain ticks to signal range (+/- a half channel)
    aud_ticks = [0.1 0.25 0.5 1 2 4 8 16 32];
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

view(90,90)
set(gca,'ygrid','on','xdir','reverse')

if save_plot
    filename = PLOT.filename;
    %             print(gcf,'-dpng',filename)
    print(gcf,'-depsc',filename)
    disp(['Plot saved as ' filename])
end

end
%EOF