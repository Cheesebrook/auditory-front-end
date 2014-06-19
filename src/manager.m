classdef manager < handle
    
    properties
        Processors      % Array of processor objects
        InputList       % Array of input "adresses" to each processors
        OutputList      % Array of output "adresses" to each processors
        Map             % Vector mapping the processing order to the 
                        % processors order. Allows for avoiding to reorder
                        % the processors array when new processors are
                        % added.
        Data            % Pointer to the data object
    end
    
    methods
        function mObj = manager(data,request,p)
            %manager        Constructs a manager object
            %
            %USAGE
            %       mObj = manager(data,request)
            %       mObj = manager(data,request,p)
            %
            %INPUT ARGUMENTS
            %     data : Handle of an existing data structure
            %  request : Cell array of requested signals, cues or features
            %            e.g., request = {'azimuth','rms'}
            %        p : Set of parameters
            %
            %OUTPUT ARGUMENTS
            %     mObj : Manager instance
            
            % TO DO:
            %  - Expand this h1 line as the implementation progresses
            %  - Add support for multiple requests

            
            if nargin>0     % Failproof for Matlab empty calls
            
            % Input check
            if nargin<3||isempty(p);p=[];end
            if nargin<2
                request = [];
            end
            if nargin<1
                error(['Too few arguments, the manager is built upon '...
                    'an existing data Object'])
            end
            
            
            % Add pointer to the data structure
            mObj.Data = data;
            
            % Instantiate the requested processors
            % TO DO: HERE GOES SUPPORT FOR MULTIPLE REQUESTS
            if ~isempty(request)
                mObj.addProcessor(request,p);
            end
            
            end
        end
        
        function processSignal(mObj)
            %processSignal      Requests a manager object to extract its
            %                   required features for a full signal present
            %                   in mObj.Data.signal
            %
            %USAGE
            %    mObj.processSignal()
            %
            %INPUT ARGUMENT
            %   mObj : Manager object
            
            % Check that there is an available signal
            if isempty(mObj.Data.signal)
                warning('No signal available for processing')
            else            
                % Number of processors
                n_proc = size(mObj.Processors,1);

                % Loop on each processor
                for ii = 1:n_proc
                    % Get index of current processor
                    jj = mObj.Map(ii);

                    if ~isprop(mObj.Processors{jj,1},'isBinaural')
                        % Apply processing for left channel (or mono if
                        % interaural cue/feature)
                        mObj.OutputList{jj,1}.Data = ...
                            mObj.Processors{jj,1}.processChunk(mObj.InputList{jj,1}.Data);

                        % Apply for right channel if stereo cue/feature
                        if mObj.Data.isStereo && ~isempty(mObj.Processors{jj,2})
                            mObj.OutputList{jj,2}.Data = ...
                                mObj.Processors{jj,2}.processChunk(mObj.InputList{jj,2}.Data);
                        end
                    else
                        % If the processor extracts a binaural cue, inputs
                        % from left and right channel should be routed
                        mObj.OutputList{jj,1}.Data = ...
                            mObj.Processors{jj,1}.processChunk(mObj.InputList{jj,1}.Data,...
                                mObj.InputList{jj,2}.Data);

                    end
                end
            end
        end
        
        function processChunk(mObj,sig_chunk)
            %processChunk   Update the signal with a new chunk of data and
            %               calls the processing chain for this new chunk
            %
            %USAGE
            %   mObj.processChunk(sig_chunk)
            %
            %INPUT ARGUMENTS
            %      mObj : Manager object
            % sig_chunk : New signal chunk
            
            % Check that the signal chunk has correct number of channels
            if size(sig_chunk,2) ~= mObj.Data.isStereo+1
                % TO DO: Change that to a warning and handle appropriately
                error(['The dimensionality of the provided signal chunk'...
                    'is incompatible with previous chunks'])
            end
            
            % Append the signal chunk
            if mObj.Data.isStereo
               mObj.Data.signal{1}.appendChunk(sig_chunk(:,1));
               mObj.Data.signal{2}.appendChunk(sig_chunk(:,2));
            else            
               mObj.Data.signal{1}.appendChunk(sig_chunk);
            end
            
            % Number of processors
            n_proc = size(mObj.Processors,1);
            
            % Loop on each processor
            for ii = 1:n_proc
                % Get index of current processor
                jj = mObj.Map(ii);
                
                if ~isprop(mObj.Processors{jj,1},'isBinaural')
                    % Apply processing for left channel (or mono if
                    % interaural cue/feature):

                    % Getting input signal handle (for code readability)
                    in = mObj.InputList{jj,1};

                    % Indexes for last chunk position in input
                    s = in.LastChunk(1);
                    e = in.LastChunk(2);

                    % Perform the processing
                    out = mObj.Processors{jj,1}.processChunk(in.Data(s:e,:,:));

                    % Store the result
                    mObj.OutputList{jj,1}.appendChunk(out);

                    % Apply similarly for right channel if binaural cue/feature
                    if mObj.Data.isStereo && ~isempty(mObj.Processors{jj,2})
                        in = mObj.InputList{jj,2};
                        s = in.LastChunk(1);
                        e = in.LastChunk(2);
                        out = mObj.Processors{jj,2}.processChunk(in.Data(s:e,:,:));
                        mObj.OutputList{jj,2}.appendChunk(out);
                    end
                    
                else
                    % Inputs from left AND right channels are needed at
                    % once
                    
                    % Getting input signal handles for both channels
                    in_l = mObj.InputList{jj,1};
                    in_r = mObj.InputList{jj,2};
                    
                    % Indexes for last chunk position in input
                    s_l = in_l.LastChunk(1);
                    e_l = in_l.LastChunk(2);
                    s_r = in_r.LastChunk(1);
                    e_r = in_r.LastChunk(2);
                    
                    % Perform the processing
                    out = mObj.Processors{jj,1}.processChunk(...
                        in_l.Data(s_l:e_l,:,:),...
                        in_r.Data(s_r:e_r,:,:));
                    
                    % Store the result
                    mObj.OutputList{jj,1}.appendChunk(out);
                    
                end
                
%                 % Getting input signal handle (for code readability)
%                 in = mObj.InputList{jj};
%                 
%                 % Indexes for last chunk position in input
%                 s = in.LastChunk(1);
%                 e = in.LastChunk(2);
%                 
%                 % Perform the processing
%                 out = mObj.Processors{jj}.processChunk(in.Data(s:e,:,:));
%                 
%                 % Store the result
%                 mObj.OutputList{jj}.appendChunk(out);
                
            end
        end
        
        function hProc = hasProcessor(mObj,name,p)
            %hasProcessor       Determines if a processor (including its
            %                   dependencies) already exists
            %
            %USAGE
            %  hProc = mObj.hasProcessor(name,p)
            %
            %INPUT ARGUMENTS
            %   mObj : Instance of manager object
            %   name : Name of processor
            %      p : Complete structure of parameters for that processor
            %
            %OUTPUT ARGUMENT
            %  hProc : Handle to an existing processor, if any, 0 else
            
            % TO DO: This function needs to recursively look into dependent
            % processors
            % - Will have to be adjusted when introducing features with
            % multiple dependencies
            
            % Initialize the output
            hProc = 0;
            
            % Loop over the processors to find the ones with suitable name
            for ii = 1:size(mObj.Processors,1)
                
                % Get a handle to that processor, for readability in the
                % following
                proc = mObj.Processors{ii,1};
                
                % Is the current processor one of the sought type?
                if isa(proc,name)
                    
                    % Does it have the requested parameters?
                    if proc.hasParameters(p)
                        
                        % Then it is a suitable candidate, we should
                        % investigate its dependencies
                        while true
                            
                            if isempty(proc.Dependencies{1})
                                % Then we reached the end of the dependency
                                % list without finding a mismatch in
                                % parameters. The original processor is a
                                % solution:
                                hProc = mObj.Processors{ii,1};
                                return
                            end
                            
                            % Set current processor to proc dependency
                            proc = proc.Dependencies{1};
                            
                            % Does the dependency also have requested
                            % parameters? If not, break of the while loop
                            if ~proc.hasParameters(p)
                                break
                            end
                            
                        end
                        
                        
                    end
                    
                end
                
                % If not, move along in the loop
                
            end
            
        end
        
        function out = addProcessor(mObj,request,p)
            %addProcessor       Add new processor needed to compute a
            %                   single request. Optionally returns a handle
            %                   for the requested signal for convenience.
            %
            %USAGE:
            %     mObj.addProcessor(request,p)
            %     out = mObj.addProcessor(...)
            %
            %INPUT ARGUMENTS
            %    mObj : Manager instance
            % request : Requested signal (string)
            %       p : Structure of non-default parameters
            %
            %OUTPUT ARGUMENTS
            %     out : Handle for the requested signal
            %
            % TO DO:
            %   - Add support for multiple requests
            %   - Add support for feedback (ie, no overwrite of existing
            %   processors.
            %   - Current bug in .cfHz property of signals. This cannot be
            %   taken from p.cfHz but needs to be fetch from dependent
            %   processors. Only affects labeling of the channels.
            
            if nargin<3 || isempty(p)
                % Initialize parameter structure
                p = struct;
            end
            
            if ~isfield(p,'fs')
                % Add sampling frequency to the parameter structure
                p.fs = mObj.Data.signal{1}.FsHz;
            end
            
            % Add default values for parameters not explicitly defined in p
            p = parseParameters(p);
            
            % Try/Catch to check that the request is valid
            try 
                % TO DO: implement for multiple requests
                getDependencies(request);
            catch err
                % Buid a list of available signals for display
                list = getDependencies('available');
                str = [];
                for ii = 1:size(list,2)-1
                    str = [str list{ii} ', '];
                end
                % Return the list
                error(['One of the requested signal, cue, or feature '...
                    'name is unknown. Valid names are as follows: %s'],str)
            end
            
            % Get the full dependency list
            if ~strcmp(request,'time')
                dep_list = [request getDependencies(request)];
            else
                % Time is a special case and is listed as its dependency
                dep_list = getDependencies(request);
            end
            
            % The processing order is the reversed list of dependencies
            dep_list = fliplr(dep_list);
            
            % Former number of processors
            n_proc = size(mObj.Processors,1);
            
            % Number of new processors involved
            n_new_proc = size(dep_list,2);
            
            % Preallocation
            if isempty(mObj.Processors)
                if mObj.Data.isStereo
                    n_chan = 2;
                else
                    n_chan = 1;
                end
                mObj.Processors = cell(n_new_proc,n_chan);   
                mObj.InputList = cell(n_new_proc,n_chan);    % TO DO: Will have to be changed to account for multiple input features
                mObj.OutputList = cell(n_new_proc,n_chan);
            end
            
            
            % Initialize pointer to dependency (first processor is always
            % timeProc)
            % TO DO: Update this to a feedback scenario, first processor is
            % the new processor with lowest dependency
            if mObj.Data.isStereo
                dep_sig_l = mObj.Data.signal{1};
                dep_sig_r = mObj.Data.signal{2};
                dep_proc_l = [];
                dep_proc_r = [];
            else
                dep_sig = mObj.Data.signal{1};
                dep_proc = [];
            end
            
            % Processors instantiation and data object property population
            for ii = n_proc+1:n_proc+n_new_proc     
                switch dep_list{ii-n_proc}
                    
                    case 'time'
                        % TO DO: Include actual time processor
                        if mObj.Data.isStereo
                            % Instantiate left and right ear processors
                            mObj.Processors{ii,1} = identityProc(p.fs);
                            mObj.Processors{ii,2} = identityProc(p.fs);
                            % Generate new signals
                            sig_l = TimeDomainSignal(mObj.Processors{ii}.FsHzOut,'time','Time',[],'left');
                            sig_r = TimeDomainSignal(mObj.Processors{ii}.FsHzOut,'time','Time',[],'right');
                            % Add the signals to the data object
                            mObj.Data.addSignal(sig_l);
                            mObj.Data.addSignal(sig_r)
                        else
                            % Instantiate a processor
                            mObj.Processors{ii} = identityProc(p.fs);
                            % Generate a new signal
                            sig = TimeDomainSignal(mObj.Processors{ii}.FsHzOut,'time','Time');
                            % Add signal to the data object
                            mObj.Data.addSignal(sig);
                        end
                                     
                    case 'gammatone'
                        if mObj.Data.isStereo
                            % Instantiate left and right ear processors
                            mObj.Processors{ii,1} = gammatoneProc(p.fs,p.f_low,p.f_high,p.IRtype,p.nERBs,p.bAlign,p.n_gamma,p.bwERBs,p.durSec);
                            mObj.Processors{ii,2} = gammatoneProc(p.fs,p.f_low,p.f_high,p.IRtype,p.nERBs,p.bAlign,p.n_gamma,p.bwERBs,p.durSec);
                            % Generate new signals
                            sig_l = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'gammatone',mObj.Processors{ii}.cfHz,'Gammatone filterbank output',[],'left');
                            sig_r = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'gammatone',mObj.Processors{ii}.cfHz,'Gammatone filterbank output',[],'right');
                            % Add the signals to the data object
                            mObj.Data.addSignal(sig_l);
                            mObj.Data.addSignal(sig_r)
                        else
                            % Instantiate a processor
                            mObj.Processors{ii} = gammatoneProc(p.fs,p.f_low,p.f_high,p.IRtype,p.nERBs,p.bAlign,p.n_gamma,p.bwERBs,p.durSec);
                            % Generate a new signal
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'gammatone',mObj.Processors{ii}.cfHz,'Gammatone filterbank output',[],'mono');
                            % Add signal to the data object
                            mObj.Data.addSignal(sig);
                        end
                        
                    case 'innerhaircell'
                        if mObj.Data.isStereo
                            % Instantiate left and right ear processors
                            mObj.Processors{ii,1} = IHCenvelopeProc(p.fs,p.IHCMethod);
                            mObj.Processors{ii,2} = IHCenvelopeProc(p.fs,p.IHCMethod);
                            % Generate new signals
                            sig_l = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'innerhaircell',p.cfHz,'Inner hair-cell envelope',[],'left');
                            sig_r = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'innerhaircell',p.cfHz,'Inner hair-cell envelope',[],'right');
                            % Add the signals to the data object
                            mObj.Data.addSignal(sig_l);
                            mObj.Data.addSignal(sig_r)
                        else
                            % Instantiate a processor
                            mObj.Processors{ii} = IHCenvelopeProc(p.fs,p.IHCMethod);
                            % Generate a new signal
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'innerhaircell',p.cfHz,'Inner hair-cell envelope',[],'mono');
                            % Add signal to the data object
                            mObj.Data.addSignal(sig);
                        end
                        
                    case 'autocorrelation'
                        if mObj.Data.isStereo
                            % Instantiate left and right ear processors
                            mObj.Processors{ii,1} = autocorrelationProc(p.fs,p);
                            mObj.Processors{ii,2} = autocorrelationProc(p.fs,p);
                            % Generate new signals
                            lags = 0:1/p.fs:mObj.Processors{ii,1}.wSizeSec-1/p.fs;   % Vector of lags
                            sig_l = CorrelationSignal(mObj.Processors{ii}.FsHzOut,'autocorrelation',p.cfHz,lags,'Auto-correlation',[],'left');
                            sig_r = CorrelationSignal(mObj.Processors{ii}.FsHzOut,'autocorrelation',p.cfHz,lags,'Auto-correlation',[],'right');
                            % Add the signals to the data object
                            mObj.Data.addSignal(sig_l);
                            mObj.Data.addSignal(sig_r)
                        else
                            % Instantiate a processor
                            mObj.Processors{ii} = autocorrelationProc(p.fs,p);
                            % Generate a new signal
                            lags = 0:1/p.fs:mObj.Processors{ii,1}.wSizeSec-1/p.fs;   % Vector of lags
                            sig = CorrelationSignal(mObj.Processors{ii}.FsHzOut,'autocorrelation',p.cfHz,lags,'Auto-correlation',[],'mono');
                            % Add signal to the data object
                            mObj.Data.addSignal(sig);
                        end
                        clear lags
                        
                    case 'crosscorrelation'
                        % Check that two channels are available
                        if ~mObj.Data.isStereo
                            warning('Manager cannot instantiate a binaural cue extractor for a single-channel signal')
                        else
                            mObj.Processors{ii} = crosscorrelationProc(p.fs,p);
                            maxLag = ceil(mObj.Processors{ii}.maxDelaySec*p.fs);
                            lags = (-maxLag:maxLag)/p.fs;
                            sig = CorrelationSignal(mObj.Processors{ii}.FsHzOut,'crosscorrelation',p.cfHz,lags,'Cross-correlation',[],'mono');
                            mObj.Data.addSignal(sig);
                            clear maxLag lags
                        end
                        
                    case 'ratemap_magnitude'
                        if mObj.Data.isStereo
                            % Instantiate left and right ear processors
                            mObj.Processors{ii,1} = ratemapProc(p.fs,p,'magnitude');
                            mObj.Processors{ii,2} = ratemapProc(p.fs,p,'magnitude');
                            % Generate new signals
                            sig_l = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ratemap_magnitude',p.cfHz,'Ratemap (magnitude)',[],'left');
                            sig_r = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ratemap_magnitude',p.cfHz,'Ratemap (magnitude)',[],'right');
                            % Add the signals to the data object
                            mObj.Data.addSignal(sig_l);
                            mObj.Data.addSignal(sig_r)
                        else
                            % Instantiate a processor
                            mObj.Processors{ii} = ratemapProc(p.fs,p,'magnitude');
                            % Generate a new signal
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ratemap_magnitude',p.cfHz,'Ratemap (magnitude)',[],'mono');
                            % Add signal to the data object
                            mObj.Data.addSignal(sig);
                        end
                        
                    case 'ratemap_power'
                        if mObj.Data.isStereo
                            % Instantiate left and right ear processors
                            mObj.Processors{ii,1} = ratemapProc(p.fs,p,'power');
                            mObj.Processors{ii,2} = ratemapProc(p.fs,p,'power');
                            % Generate new signals
                            sig_l = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ratemap_power',p.cfHz,'Ratemap (power)',[],'left');
                            sig_r = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ratemap_power',p.cfHz,'Ratemap (power)',[],'right');
                            % Add the signals to the data object
                            mObj.Data.addSignal(sig_l);
                            mObj.Data.addSignal(sig_r)
                        else
                            % Instantiate a processor
                            mObj.Processors{ii} = ratemapProc(p.fs,p,'power');
                            % Generate a new signal
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ratemap_power',p.cfHz,'Ratemap (power)',[],'mono');
                            % Add signal to the data object
                            mObj.Data.addSignal(sig);
                        end
                        
                    case 'ild'
                        % Check that two channels are available
                        if ~mObj.Data.isStereo
                            warning('Manager cannot instantiate a binaural cue extractor for a single-channel signal')
                        else
                            mObj.Processors{ii} = ildProc(p.fs,p);
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ild',p.cfHz,'Interaural Level Difference',[],'mono');
                            mObj.Data.addSignal(sig);
                        end
                        
                    case 'ic_xcorr'
                        if ~mObj.Data.isStereo
                            warning('Manager cannot instantiate a binaural cue extractor for a single-channel signal')
                        else
                            mObj.Processors{ii} = icProc(dep_proc_l.FsHzOut,p);
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'ic_xcorr',p.cfHz,'Interaural correlation',[],'mono');
                            mObj.Data.addSignal(sig);
                        end
                        
                    case 'itd_xcorr'
                        if ~mObj.Data.isStereo
                            warning('Manager cannot instantiate a binaural cue extractor for a single-channel signal')
                        else
                            mObj.Processors{ii} = itdProc(dep_proc_l.FsHzOut,p);
                            sig = TimeFrequencySignal(mObj.Processors{ii}.FsHzOut,'itd_xcorr',p.cfHz,'Interaural Time Difference',[],'mono');
                            mObj.Data.addSignal(sig);
                        end
                        
                    % TO DO: Populate that list further
                    
                    % N.B: No need for "otherwise" case once complete
                    
                    otherwise
                        error('%s is not supported at the moment',...
                            dep_list{ii+1});
                end
                
                
                % Add input/output pointers, dependencies, and update dependencies.
                % Three possible scenarios:
                
                if isprop(mObj.Processors{ii},'isBinaural')
                    
                    % 1-Then there are two inputs (left&right) and one output
                    mObj.InputList{ii,1} = dep_sig_l;
                    mObj.InputList{ii,2} = dep_sig_r;
                    mObj.OutputList{ii,1} = sig;
                    mObj.OutputList{ii,2} = [];
                    mObj.Processors{ii,1}.Dependencies = {dep_proc_l,dep_proc_r};
                    dep_sig = sig;
                    dep_proc = mObj.Processors{ii};
                    
                elseif exist('sig','var')&&strcmp(sig.Canal,'mono')
                    
                    % 2-Then there is a single input and single output
                    mObj.InputList{ii,1} = dep_sig;
                    mObj.OutputList{ii,1} = sig;
                    mObj.Processors{ii}.Dependencies = {dep_proc};
                    dep_sig = sig;
                    dep_proc = mObj.Processors{ii};
                    
                else
                    
                    % 3-Else there are two inputs and two outputs
                    mObj.InputList{ii,1} = dep_sig_l;
                    mObj.InputList{ii,2} = dep_sig_r;
                    mObj.OutputList{ii,1} = sig_l;
                    mObj.OutputList{ii,2} = sig_r;
                    mObj.Processors{ii,1}.Dependencies = {dep_proc_l};
                    mObj.Processors{ii,2}.Dependencies = {dep_proc_r};
                    dep_sig_l = sig_l;
                    dep_sig_r = sig_r;
                    dep_proc_l = mObj.Processors{ii,1};
                    dep_proc_r = mObj.Processors{ii,2};
                    
                end
                

                % Clear temporary handles to ensure no inconsistencies 
                clear sig sig_l sig_r
                
            end
            
            % The mapping at this point is linear
            mObj.Map(n_proc+1:n_proc+n_new_proc) = n_proc+1:n_proc+n_new_proc;
            
            % Provide the user with a pointer to the requested signal
            if nargout>0
                out = mObj.OutputList{n_proc+n_new_proc,1};
            end
            
        end
        
        function [hProc,list] = findInitProc(mObj,request,p)
            %findInitProc   Find an initial compatible processor for a new
            %               request
            %
            %USAGE:
            %         hProc = mObj.findInitProc(request,p)
            %  [hProc,list] = mObj.findInitProc(request,p)
            %
            %INPUT PARAMETERS
            %    mObj : Manager instance
            % request : Requested signal name
            %       p : Parameter structure associated to the request
            %
            %OUTPUT PARAMETERS
            %   hProc : Handle to the highest processor in the processing 
            %           chain that is compatible with the provided
            %           parameters
            %    list : List of signal names that need to be computed,
            %           starting from the output of hProc, to obtain the
            %           request
        
            % Input parameter checking
            if nargin<3 || isempty(p)
                % Initialize parameter structure
                p = struct;
            end
            if ~isfield(p,'fs')
                % Add sampling frequency to the parameter structure
                p.fs = mObj.Data.signal{1}.FsHz;
            end
            % Add default values for parameters not explicitly defined in p
            p = parseParameters(p);
        
            % Try/Catch to check that the request is valid
            try
                getDependencies(request);
            catch err
                % Buid a list of available signals for display
                list = getDependencies('available');
                str = [];
                for ii = 1:size(list,2)-1
                    str = [str list{ii} ', '];
                end
                % Return the list
                error(['The requested signal, %s is unknown. '...
                    'Valid names are as follows: %s'],request,str)
            end
            
            % Get the full list of dependencies corresponding to the request
            if ~strcmp(request,'time')
                dep_list = [request getDependencies(request)];
            else
                % Time is a special case as it is listed as its own dependency
                dep_list = getDependencies(request);
            end
            
            % Initialization of while loop
            ii = 1;
            dep = signal2procName(dep_list{ii});
            hProc = mObj.hasProcessor(dep,p);
            list = {};
            
            % Looping until we find a suitable processor in the list of
            % dependency
            while hProc == 0 && ii<size(dep_list,2)
                
                % Then we will need to re-compute that signal
                list = [list dep_list{ii}];
                
                % Move on to next level of dependency
                ii = ii + 1;
                dep = signal2procName(dep_list{ii});
                hProc = mObj.hasProcessor(dep,p);
                
            end
            
            if hProc == 0
                % Then all the signals need recomputation, including time
                list = [list dep_list{end}];
            end
            
        end
        
    end
    
    
end