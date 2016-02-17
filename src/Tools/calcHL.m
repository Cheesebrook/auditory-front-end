%calcHL.m
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
% This processor estimates the Hearing Level (HL) either based on 
% the absolute threshold of hearing  ('ath') or on the energy in each filter ('en').
% The thresholds are computed for each filter of the auditory filterbank.
% Implementation as in hlProc.m!
%
% inputs:
% in        - Left and right input signal [time x frequency]
% method    - Method used to estimate the thresholds 
%             absolute threshold of hearing ('ath') or energy-based ('en')
% enTh      - Threshold for energy-based method
%
% output:
% out       - Left and right input signal with applied threshold
% totEn     - SPL in each frequency-band ('ath')
%             Total energy in each frequency-band ('en')
% actbands  - Vector containing the frequency-bands that are above
%             threshold
%
% --------History-----------
% - 30. September 2015: Modelling ASW framework
% - 23. October 2015
% --------------------------

function [out_l, out_r, totEn, binMask, freqWeights] = ...,
         calcHL(in_l,in_r,method,cfHz,FsHzIn,bframebased,enTh,hl_wname,hl_wSizeSec,hl_hSizeSec)

%% Check Inputs
if nargin<2
    error('Please provide left and right input!')
end
if nargin<3 || isempty(method)
    method = 'none';
end
if nargin<4 || isempty(FsHzIn)
    error('Please provide a center-frequency vector!')
end
if nargin<5 || isempty(FsHzIn)
    FsHzIn = 44.1e3; %Sampling Frequency [Hz]
end
if nargin<6|| isempty(bframebased)
    bframebased = 1; %bool: frame-based estimation(1)
end
if nargin<7 || isempty(enTh)
    enTh = 0.5; %Threshold for energy-estimation
end
if nargin<8 || isempty(hl_wname)
    hl_wname = 'hann'; %Window name
end
if nargin<9|| isempty(hl_wSizeSec)
    hl_wSizeSec = 20E-3; %Window duration (s)
end
if nargin<10|| isempty(hl_hSizeSec)
    hl_hSizeSec = 10E-3; %Window step size (s)
end

if ~bframebased&&strcmp(method,'binmask')
    error(['The binary mask method can obiously only be applied ',...
        'for frame-based processing! Set bframebased to ''true''!'])
end

% get dimensions of input
[nSamples,nChannels] = size(in_l);

% Window
wSize = 2*round(hl_wSizeSec*FsHzIn/2);
hSize = round(hl_hSizeSec*FsHzIn);
win = window(hl_wname,wSize);

% How many frames are in the buffered input?
nFrames = floor((nSamples-(wSize-hSize))/hSize);

% Compute HL:

% Pre-allocate outputs
out_l = NaN(nSamples,nChannels);
out_r = NaN(nSamples,nChannels);

% init SPL calculations
frame_l_rms = zeros(nFrames,nChannels);
frame_r_rms = zeros(nFrames,nChannels);

% Window processing
winRep = repmat(win,1,nChannels);

% Estimate SPL
xref = 20e-6; %ref for dB SPL definition: xrms = 1 -> Lp = 94dB SPL
%(xrms = xref = 20e-6 -> Lp = 0dB SPL)
% xref = 1e-5; %ref for digital value definition: xrms = 1 -> Lp = 100dB
%              %using this reference will reproduce totEnergy = Lp
%              (Lp being the specified level of the input signal)

if bframebased
    % Loop over time frames
    for ii = 1:nFrames
        % Get start and end indexes for the current frame
        n_start = (ii-1)*hSize+1;
        n_end = (ii-1)*hSize+wSize;
        
        % Energy in the windowed frame for left and right input
        frame_l_rms(ii,:) = sqrt(mean(power(winRep.*in_l(n_start:n_end,:),2)));
        frame_r_rms(ii,:) = sqrt(mean(power(winRep.*in_r(n_start:n_end,:),2)));
    end
    
    % calc SPL and devide by number of frames for averaging
    frame_l_SPL = 20*log10(1/nFrames*sum(frame_l_rms/xref,1)); %acc. SPL left channel
    frame_r_SPL = 20*log10(1/nFrames*sum(frame_r_rms/xref,1)); %acc. SPL right channel
    
else
    % Energy in the windowed frame for left and right input
    frame_l_rms = sqrt(mean(power(in_l,2)));
    frame_r_rms = sqrt(mean(power(in_r,2)));
    
    frame_l_SPL = 20*log10(frame_l_rms/xref); %acc. SPL left channel
    frame_r_SPL = 20*log10(frame_r_rms/xref); %acc. SPL right channel
end

% select maximum values from left and right channel
max_lr_SPL = max(frame_l_SPL,frame_r_SPL);

% frame_SPL contains total SPL after frame-based processing
% (check if it corresponds to the specified level)
totEn = 20*log10(sum(10.^(max_lr_SPL/20))); %tot. level from all filters

% Weights normalized to frequency band with highest SPL
freqWeights = max_lr_SPL/max(max_lr_SPL);
freqWeights(freqWeights<0) = 0; %set negative weights (due to subthreshold) to zero

% calculate active frequency bands
switch method
    case 'none'
        %binary mask: all time-frequency-units active
        binMask = ones(nFrames,nChannels);
        
    case 'ath'
        %absolute threshold of hearing
        ath = calcATH(cfHz*1e-3); 

        %all filters that are above absolute threshold of hearing
        bpassth = max_lr_SPL>ath; %bool
        
        %binary mask: all time-frequency-units active
        binMask = repmat(bpassth,nFrames,1);
        
    case 'binmask'
        %absolute threshold of hearing
        ath = calcATH(cfHz*1e-3); 
        
        % maximal time-frequency-units between left and right channel
        frame_SPL = max(20*log10(frame_l_rms/xref),20*log10(frame_r_rms/xref));
        
        %all filters that are above absolute threshold of hearing
        binMask = frame_SPL>repmat(ath,nFrames,1); %time-frequency binary mask 
        
%         %binary mask: "substitutes" actbands
%         binMask = ones(size(frame_SPL)); %column-vector with #filters
%         binMask = binMask(bpassth); %select activated channels from the vector
        
        % dynamic time-frequency weights
        % timeFreqWeights = frame_SPL/max(max(frame_SPL));
        % timeFreqWeights(timeFreqWeights<0) = 0; %set negative weights (due to subthreshold) to zero
        
        case 'binmask30'
        %absolute threshold of hearing
        ath = calcATH(cfHz*1e-3);
        offset = 30; 
        
        % maximal time-frequency-units between left and right channel
        frame_SPL = max(20*log10(frame_l_rms/xref),20*log10(frame_r_rms/xref));
        
        %all filters that are 30 dB above absolute threshold of hearing
        binMask = frame_SPL>repmat(ath+offset,nFrames,1); %time-frequency binary mask 
        
    case 'en'
        warning('Energy-based method not implemented yet!')
        
    otherwise
        error(['The method ' pObj.Method ' is not defined!'])
end

% Set all samples to 0 that are below the specified threshold
for ii = 1:nFrames
    % Get start and end indexes for the current frame
    n_start = (ii-1)*hSize+1;
    n_end = (ii-1)*hSize+wSize;
    ncurrSamples = n_end-n_start+1;
    
    out_l(n_start:n_end,:) = winRep.*in_l(n_start:n_end,:).*repmat(binMask(ii,:),ncurrSamples,1);
    out_r(n_start:n_end,:) = winRep.*in_r(n_start:n_end,:).*repmat(binMask(ii,:),ncurrSamples,1);
end

% Histogram of thresholded SPL across all frames
%hist(frame_SPL(binMask))
% hist(frame_SPL(:))
% xlabel('SPL'), ylabel('Counts'), title('DS2014')

%EOF