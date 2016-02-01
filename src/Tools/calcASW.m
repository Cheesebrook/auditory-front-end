%calcASW.m
%
% asw = (dObj,varargin)
%
% This function predicts apparent source width (ASW) perception from ITD and ILDs
% !!!In future this function should be replaced by a the corresponding
% processor aswProc.m
%
% Inputs:
% dObj - data object containing ITDs and ILDs
% 
% optional inputs (varargin):
% wdMethod - Method for estimating the distribution width ('prct' or 'std')
% percent  - Choose percent values for prctile calculation
% avMethod - Method for averaging across frequency channels
% transform - 'Choose transformation of bin. cue: 'none', 'norm' (normalize)  or 'latcomp' (lateral compression)'
% freqWeighting - Select frequency-dependent weighting of the binaural cue
%
% Outputs:
% asw - Prediction of ASW perception.
%
% by Johanens Käsbach, CAHR, DTU, johk@elektro.dtu.dk, 10. September 2015
%
%---------History-----------
% - no further editing done
%
%---------------------------

function [asw, aswReprChan, itdReprChan, ildReprChan] = calcASW(dObj,varargin)

%% Inputs
% check inputs
if nargin<1||isempty(dObj)
    error('Please provide a data object!')
end
% if isempty(findprop(dObj,'itd'))||isempty(findprop(dObj,'ild'))
%     error('Please provide an object containing property ''itd'' and ''ild''!')
% end
% temporarily use this implementation (due to substituted hlProc.m)
if ~isfield(dObj,'itd')||~isfield(dObj,'ild')
    error('Please provide an object containing property ''itd'' and ''ild''!')
end

% optional parameters
opt = struct(varargin{:});
                             
if isfield(opt,'wdMethod')
    wdMethod = opt.wdMethod;
else
    wdMethod = 'prct';
end
if isfield(opt,'percent')
    percent = opt.percent;
else
    percent = [10 90];
end
if isfield(opt,'itdmax')
    itdmax = opt.itdmax;
else
    itdmax = 1e-3;
end
if isfield(opt,'ildmax')
    ildmax = opt.ildmax;
else
    ildmax = 20;
end
if isfield(opt,'chanRepresent')
    chanRepresent = opt.chanRepresent;
else
    chanRepresent = 'width';
end
if isfield(opt,'transMethod')
    transMethod = opt.transMethod;
else
    transMethod = 'none';
end
if isfield(opt,'freqWeighting')
    freqWeighting = opt.freqWeighting;
else
    freqWeighting = 'none';
end
if isfield(opt,'combMethod')
    combMethod = opt.combMethod;
else
    combMethod = 'itd';
end
if isfield(opt,'freqWeights')
    freqWeights = opt.freqWeights;
else
    freqWeights = ones(1,size(dObj.itd,2)); %use ones, so no weighting!
end

%% Definitions

%% Processing
% get itd and ild data
% itd = dObj.itd{1}.Data(:);
% ild = dObj.ild{1}.Data(:);
% cfHz = dObj.itd{1}.cfHz;
% Nchan = size(dObj.itd{1}.Data(:),2);
% temporarily use this implementation (due to substituted hlProc.m)
itd = dObj.itd;
ild = dObj.ild;
cfHz = dObj.cfHz;
Nchan = size(dObj.itd,2);

% Width per channel
[itdWidthChan, itdLR_boarderChan, itdPrctChan] = calcDistrWidth(itd,wdMethod,percent); %calculate width of binCue directly for each frequency channel
[ildWidthChan, ildLR_boarderChan, ildPrctChan] = calcDistrWidth(ild,wdMethod,percent); %calculate width of binCue directly for each frequency channel

% choose representation of ASW per channel
switch chanRepresent
    case 'width' %width
        itdReprChan = itdWidthChan;
        ildReprChan = ildWidthChan;
        
    case 'lr' %left and right boundary seperately
        itdReprChan = itdLR_boarderChan;
        ildReprChan = ildLR_boarderChan;
        
    case 'prct' %use all percentiles
        if strcmp(chanRepresent,wdMethod)
            itdReprChan = itdPrctChan;
            ildReprChan = ildPrctChan;
        else
            warning(['This representation method is only valid for a '],...
                  ['percentile-based calculation! Using ''lr'' method instead!'])
            itdReprChan = itdLR_boarderChan;
            ildReprChan = ildLR_boarderChan;             
        end  
end

% Transformation of bin. cue
switch transMethod
    case 'none' %no transformation
        %do nothing
        
    case 'norm' %normalization
        itdReprChan = itdReprChan/itdmax;
        ildReprChan = ildReprChan/ildmax;
        
    case 'freqnorm'
        % normalize to the highest possible value from anechoic 
        % conditions in each band
        itdload  = load('ITD2Azimuth_Subband.mat');
        itdmax   = max(abs(itdload.mapping.itd2azim),[],1);
        itdReprChan = itdReprChan.*repmat(itdmax.^-1,[size(itdReprChan,1),1]);
        ildload  = load('ILD2Azimuth_Subband.mat');
        ildmax   = max(abs(ildload.mapping.ild2azim),[],1);
        ildReprChan = ildReprChan.*repmat(ildmax.^-1,[size(itdReprChan,1),1]);
        
    case 'latcomp' %lateral compression
        itdReprChan = latCompression(itdReprChan,itdmax);
        ildReprChan = latCompression(ildReprChan,ildmax);
        
    case 'mapping'  %Map boarders (percentiles/std) of ITDs and ILDs to azimuthal angle
        % init
        Nbounds = size(itdReprChan(:,1),1);
        itd_mapped = zeros(Nbounds,Nchan);
        ild_mapped = zeros(Nbounds,Nchan);
        % angleoffset = -91; %offset to find the correct angle [degree]
        itdload = load('ITD2Azimuth_Subband.matdau');
        itd2azim = itdload.mapping.itd2azim;
        ildload = load('ILD2Azimuth_Subband.mat');
        ild2azim = ildload.mapping.ild2azim;
        
        for ii = 1:Nchan %all channels
            % itd mapping
            azimuth = -90:90;
            itd_mapped(:,ii) = interp1(itd2azim(:,ii),azimuth,itdReprChan(:,ii),'nearest','extrap');
            
            % ild mapping
            ild_mapped(:,ii) = interp1(ild2azim(:,ii),azimuth,ildReprChan(:,ii),'nearest','extrap');
        end
        
        % Overwrite channel representation per cue with results from
        % mapping
        itdReprChan = itd_mapped;
        ildReprChan = ild_mapped;

otherwise
    error(['The transformation method ' transMethod ' is not defined!'])
    
end

% Frequency weighting
switch freqWeighting
    case 'none' %no weighting, i.e. use all bands 
        itdBandSelect = 1:Nchan;
        ildBandSelect = 1:Nchan;
        
    case 'itdlow' %lowpass itds, i.e. use itds only at low freqeuncies
        itdBandSelect = 1:19; %lowpass filter cut-off [#frequency band]
        ildBandSelect = 1:Nchan;
        
    case 'itdE3' %lowpass itds, i.e. use itds only at low freqeuncies
        itdBandSelect = 7:22; %#frequency bands that correspond to the averaging 
                        %performed in IACC_E3, 
                        %i.e. octave bands at    0.5,   1   and 2 kHz
                        %corresponding to bands 7:12, 11:16 and 17:22,
                        %respectively
        ildBandSelect = 1:Nchan;
        
    case 'spl'
        itdBandSelect = 1:Nchan;
        ildBandSelect = 1:Nchan;
        itdReprChan = itdReprChan.*repmat(freqWeights,size(itdReprChan,1),1);
        ildReprChan = ildReprChan.*repmat(freqWeights,size(ildReprChan,1),1);
        
%     case 'timefreqspl'
    case 'bands';
        % form matrix with 39 single bands plus 39 accumulated band steps
        % treat NAN bands
        Nbands = size(itdReprChan,2);
        
%         BandSelect = eye(2*Nbands);
%         % looking for the most significant: 18,19 plus x
%         BandSelect([17 18 19],:) = 1; %
%         % norm
%         BandSelect = [eye(2*Nbands) BandSelect * diag(1./sum(BandSelect,1))];
%         imagesc(BandSelect)
        
        itdBandSelect = [eye(Nbands) triu(ones(Nbands)')];
        ildBandSelect = itdBandSelect;
        imagesc(itdBandSelect)

        % prevent nan results
        nanBands = find(isnan(itdReprChan(1,:)));
        itdReprChan(:,nanBands) = 0;
        ildReprChan(:,nanBands) = 0;
        
    otherwise
        error(['The frequency weighting ' freqWeighting ' is not defined!'])
end

% Combination of both cues
switch combMethod
    case 'itd' %use only itd
        aswReprChan = itdReprChan(:,itdBandSelect);
        asw = nanmean(itdReprChan(:,itdBandSelect),2); %average all channels
        
    case 'ild' %use only ild
        aswReprChan = ildReprChan(:,ildBandSelect);
        asw = nanmean(ildReprChan(:,ildBandSelect),2); %average all channels
        
    case 'duplex' %combine itd and ild according to duplex theory
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be combined!')
        end
        
        % cross-over frequency in Hz
        fcrossHz = 1500;
        % find corresponding channel
        fcrossBand = interp1(cfHz,1:Nchan,fcrossHz,'next','extrap');
        
        % init
        aswReprChan = zeros(size(itdReprChan));
        
        % choose channels below fcrossHz from itds and above from ilds
        aswReprChan(:,1:fcrossBand) = itdReprChan(:,1:fcrossBand);
        aswReprChan(:,fcrossBand+1:end) = ildReprChan(:,fcrossBand+1:end);
         
%         aswReprChan(:,[3]) = itdReprChan(:,[3]);
%         aswReprChan(:,[18]) = ildReprChan(:,[18]);
%         aswReprChan(aswReprChan == 0) = nan;
        
        % average all channels
        asw = nanmean(aswReprChan,2);
        
    case 'dominant' %choose the dominant cue in each channel, i.e. the one with higher variance
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be combined!')
        end
        
        bitdDominance = (itdReprChan>ildReprChan);
        
        % init
        aswReprChan = zeros(size(itdReprChan));
        
        % choose channels for either dominance
        aswReprChan(bitdDominance) = itdReprChan(bitdDominance);
        aswReprChan(~bitdDominance) = ildReprChan(~bitdDominance);  
        
        % average all channels
        asw = nanmean(aswReprChan,2);
        
    case 'perimeter' % choose the outer bands, as done by Ahrens
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be combined!')
        end
        
        periITD = find(cfHz<200);
        periILD = find((cfHz>2e3));
        
        % init
        %aswReprChan = zeros(size(itdReprChan));
        aswReprChan = zeros(size(itdReprChan,1),size([periITD,periILD],2));
        
        % choose channels for either dominance
        aswReprChan(:,periITD)           = itdReprChan(:,periITD); %1:4
        aswReprChan(:,end-periILD+1:end) = ildReprChan(:,periILD); %5:24
        
        % average all channels
        asw = nanmean(aswReprChan,2);
    case 'perimeter2' % choose the outer bands, as done by Ahrens
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be combined!')
        end
        
        periITD = find((200<cfHz) & (cfHz<500));
        periILD = find((4e3<cfHz) & (cfHz<8e3));
        
        % init
        %aswReprChan = zeros(size(itdReprChan));
        aswReprChan = zeros(size(itdReprChan,1),size([periITD,periILD],2));
        
        % choose channels for either dominance
        aswReprChan(:,1:length(periITD))         = itdReprChan(:,periITD); %1:4
        aswReprChan(:,end-length(periILD)+1:end) = ildReprChan(:,periILD); %5:24
        
        % average all channels
        asw = nanmean(aswReprChan,2); 
        
    case 'successive'
        if ~strcmp(freqWeighting,'bands')
            warning('combMethod "successive" expects freqWeighting "bands"')
        end
         asw = [itdReprChan*itdBandSelect ildReprChan*ildBandSelect];
         aswReprChan = asw;
%         aswReprChan = [itdReprChan ildReprChan]*BandSelect;
%         asw = [itdReprChan ildReprChan]*BandSelect;
        asw(asw==0) = nan;
    otherwise
        error(['The combination method ' combMethod ' is not defined!'])

end

%EOF