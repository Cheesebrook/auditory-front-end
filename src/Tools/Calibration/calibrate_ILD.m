function calibrate_ILD(fs,winSec,cfHz,bCalibrate,wdMethod,percent)
%calibrate_ILD_Subband   Calculate frequency-dependent ILD2azimuth mapping
%
%USAGE
%      calibrate_ILD(fs,winSec,cfHz)
%      calibrate_ILD(fs,winSec,cfHz,bCalib)
%
%INPUT PARAMETERS
%            fs : sampling frequency in Hertz
%        winSec : frame size in seconds of the cross-correlation analysis
%          cfHz : Vector of the gammatone filterbank center frequencies. Empty if
%                 broadband
%    bCalibrate : if true, enforce re-computation of the mapping function
% 
%OUTPUT PARAMETERS
%     The ILD2Azimuth mapping will be stored in the MAT file
%     ILD2Azimuth_Subband.mat or ILD2Azimuth_Subband.mat
%     inside the \Tools directory.
%
%   Developed with Matlab 8.2.0.701 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/22
%   v.0.2   2014/01/28 added flag to enforce re-calibration 
%   v.0.3   2015/08/04 unified broadband and subband and adapted for using AFE (RD)
%           2016/01/22 added estimated ITD fluctuations via specified
%           wdMethod (as performed in calcASW.m) by Johannes Käsbach (JK)
%   ***********************************************************************

% Initialize persistent variables
persistent PERfs PERwinSec PERcfHz


%% 1. CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 3 || nargin > 6
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 4 || isempty(bCalibrate); bCalibrate = false; end

if nargin < 5 || isempty(wdMethod)
   wdMethod = 'prct';
end

if nargin < 6 || isempty(percent)
    percent = [10 30 50 70 90];
end

wdMethod = 'std';

%% 2. CALIBRATION SETTINGS
% 
% 
% Length of random noise sequence in seconds
lengthSec = 1;

% Use anechoic HRTF measurements
room = 'SURREY_A';
% room = 'SURREY_ROOM_B';

% Azimuth range of interest (real sound source positions)
azimRange = (-90:5:90);

% New azimuth range after interpolation
azimRangeInterp = (-90:1:90);

% Order of polynomial fit that is applied to the ILD to ensure a monotonic
% mapping
pOrder = 3;
    

%% 3. CALIBRATION STAGE
% 
% 
% Check if we can re-use the calibration file from the last function call
if isequal(fs,PERfs) && isequal(winSec,PERwinSec) && ...
   (isequal(cfHz,PERcfHz) || isempty(cfHz)) && ~bCalibrate 
    bRecalibrate = false;
elseif ~bCalibrate
    bRecalibrate = true;
    
    % Check if the existing calibration files can be used
    try     % Try/catch for particular cases (e.g., no existing calibration files
        if isempty(cfHz)
            % Broadband
            load([pwd,filesep,'src',filesep,'Summer_school_2015',filesep,'Tools',filesep,'ILD2Azimuth_Broadband.mat'],'mapping');    
            if isequal(winSec,mapping.winSec) && isequal(fs,mapping.fs) %#ok<NODEF>
                bRecalibrate = false;
            end
        else
            % Subband
            load([pwd,filesep,'src',filesep,'Summer_school_2015',filesep,'Tools',filesep,'ILD2Azimuth_Subband.mat'],'mapping');    
            if isequal(winSec,mapping.winSec) && isequal(fs,mapping.fs) ...
                    && isequal(cfHz,mapping.cfHz) %#ok<NODEF>
                bRecalibrate = false;
            end
        end
    catch
        
    end
else
    bRecalibrate = true;
end

% Perform calibration
if bRecalibrate
    % Store persistent variables
    PERfs = fs; PERwinSec = winSec;
    if ~isempty(cfHz)
        PERcfHz = cfHz;
    end

    % Number of different sound source positions
    nAzim = numel(azimRange);
    
    % Number of different sound source positions after interpolation
    nAzimInterp = numel(azimRangeInterp);
    
    % Create white noise
    noise = randn(round(lengthSec*fs),1);

    % Number of aduitory filters
    if ~isempty(cfHz)
        nFilter = numel(cfHz);
    else
        nFilter = 1;
    end
    
    % Allocate memory
    ild2Azim       = zeros(nAzim,nFilter);
    ild2AzimInterp = zeros(nAzimInterp,nFilter);
    ild2AzimPoly   = zeros(nAzimInterp,nFilter);
    ild2AzimFluct  = zeros(nAzim,nFilter);
    ild2AzimFluctInterp = zeros(nAzimInterp,nFilter);
    ild2AzimFluctPoly   = zeros(nAzimInterp,nFilter);
            
    % Initialize the AFE for subband ILD measurement
    if ~isempty(cfHz)
        par = genParStruct('cc_bBroadband',0,'cc_wSizeSec',winSec,...
                           'cc_hSizeSec',winSec/2,'cc_maxDelaySec',1.25e-3,...
                           'fb_cfHz',cfHz,'fb_nERBs',1,...
                           'ihc_method','none');
    else
        par = genParStruct('cc_bBroadband',1,'cc_wSizeSec',winSec,...
                           'cc_hSizeSec',winSec/2,'cc_maxDelaySec',1.25e-3);
    end
    dObj = dataObject([],fs,10,2);
    mObj = manager(dObj,'ild',par);
    
    if isempty(cfHz)
        name = 'Broadband';
    else
        name = 'Subband';
    end
    
    % MAIN LOOP
    %
    %
    % Loop over number of different sound source directions
    for ii = 1 : nAzim

        % Spatialize audio signal using HRTF processing
        binaural = spatializeAudio(noise,fs,azimRange(ii),room);
        
        % Estimate ILD
        mObj.processSignal(binaural)
        ildEst = mean(dObj.ild{1}.Data(:),1);
        
         % Estimate ILD fluctuation via specified wdMethod
        ildWidthChan = calcDistrWidth(dObj.ild{1}.Data(:),wdMethod,percent); %calculate width of binCue directly for each frequency channel

        % Store azimuth-dependent ILD
        ild2Azim(ii,:) = ildEst;
        
        % Store azimuth-dependent ITD fluctuation
        ild2AzimFluct(ii,:) = ildWidthChan;
        
        % Report progress
        fprintf('\n%s-based ILD2Azimuth calibration: %.2f %%',name,100*ii/nAzim);
    end
    
    fprintf('\n')
        
    % Interpolation
    %
    %
    % Loop over the number of files
    for jj = 1 : nFilter
        % Interpolate to 'rangeAzInterp'
        ild2AzimInterp(:,jj) = interp1(azimRange,ild2Azim(:,jj),azimRangeInterp);
        
        % Ensure that mapping is monotonic by using a polynomial fit
        ild2AzimPoly(:,jj) = polyval(polyfit(azimRangeInterp,ild2AzimInterp(:,jj).',pOrder),azimRangeInterp);
        
        % Interpolate to 'rangeAzInterp'
        ild2AzimFluctInterp(:,jj) = interp1(azimRange,ild2AzimFluct(:,jj),azimRangeInterp);
        
        % Ensure that mapping is monotonic by using a polynomial fit
        ild2AzimFluctPoly(:,jj) = polyval(polyfit(azimRangeInterp,ild2AzimFluctInterp(:,jj).',pOrder),azimRangeInterp);
    end
    
    
    % Save data
    %
    %
    mapping.fs           = fs;
    mapping.azim         = azimRangeInterp;
%    mapping.ild          = dObj.crosscorrelation{1}.lags;
    mapping.ild2azimRaw  = ild2Azim;
    mapping.ild2azim     = ild2AzimPoly;
    mapping.ild2azimFluctRaw = ild2AzimFluct;
    mapping.ild2azimFluct = ild2AzimFluctInterp;
    mapping.polyOrder    = pOrder;
    mapping.ildMax       = max(ild2AzimPoly);
    mapping.ildMin       = min(ild2AzimPoly); 
    mapping.cfHz         = cfHz;
    mapping.winSec       = winSec;

    % Store ILD2Azimuth template
    save([fullfile(fileparts(which(mfilename))),filesep,'ILD2Azimuth_' name '.mat'],'mapping');
end