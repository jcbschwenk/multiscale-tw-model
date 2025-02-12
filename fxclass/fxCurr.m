classdef fxCurr < matlab.mixin.Copyable
    properties
        I(1,:,:) double % dim order [1 x nTm x nTr]
        isOn(1,:,:) double
        t(1,:) double 
        Name(1,:) char
    end
    methods
        function obj = fxCurr(varargin)
             
            p = inputParser;
            p.addParameter('I', []);
            p.addParameter('Name', []);
            
            p.parse(varargin{:})
            
            obj.I = p.Results.I;
            obj.Name = p.Results.Name;
        end
        %%
        function out = isNull(obj)
            for i = 1:numel(obj)
                out(i) = isempty(obj(i).I) || ~any(obj(i).I(:));
            end
        end
        %%
        function out = dur(obj)
            out = size(obj.I,2);
        end
        %%
        function n = numTrials(obj)
            n = size(obj.I,3);
        end
        %%
        function obj = duplicateTrials(obj, nTr)
            for runObj = unique(obj)
                assert(runObj.numTrials == 1, 'Object already has > 1 trial');
                runObj.I(1,:,2:nTr) = repmat(runObj.I(1,:,1), 1,1,nTr - 1);
            end
        end
        %%
        function obj = set(obj, type, varargin)
            p = inputParser;
            
            p.addParameter('Trial', []); % which trial to set, all if empty
            
            % shared parameters:
            p.addParameter('Time', []); % pass predefined time vector (e.g. from another fxCurr obj)
            p.addParameter('Dur', 1e3); 
            p.addParameter('Pre', 200); 
            p.addParameter('Baseline', 0); % corresponding to off state
            p.addParameter('On', [-Inf Inf]); % on times
            p.addParameter('Scale', 1); % scaling of on-states, can be array for each win
            
            % for 'SmoothPulse':
            p.addParameter('RiseKernel', 20); % corresponding to sigma of convoluting gaussian
            p.addParameter('FallKernel', 20);
            
            % for 'Sine':
            p.addParameter('Frequency', 10); % time step 1ms assumed
            p.addParameter('PhaseOffset', 0); % in radians
    
            % for 'WhiteNoise':
            p.addParameter('FixedSeed', false); % if true, noise will be identical across trials
            
            if numel(varargin)==1 &&...
                    ischar(varargin{1}) &&...
                    strcmpi(varargin{1}, 'help')
                p.parse();
                disp(p.Results);
                return
            else
                p.parse(varargin{:});
            end
            
            
            
            
            whichTrial = p.Results.Trial;
            if isempty(whichTrial)
               whichTrial = 1:max(1, obj.numTrials); 
            end
            nTr = numel(whichTrial);
                        
            on = p.Results.On;
            if ~iscell(on)
                on = {on};
            end
            
            scale = p.Results.Scale;
            if isscalar(scale)
                scale = ones(1, numel(on)) .* scale;
            end
            
            if isempty(p.Results.Time)
                dur = p.Results.Dur + p.Results.Pre;
                tm = (1:dur) - p.Results.Pre;
            else
                tm = p.Results.Time;
                dur = numel(tm);
            end
            
            out = ones(1,dur).*p.Results.Baseline;
            for ievent = 1:numel(on)
                out(obj.rngfn(tm, on{ievent})) = scale(ievent) + p.Results.Baseline;
            end
            
            switch type
                
                case {'Pulse', 'Rect', 'RectPulse'}
                    
                    out = repmat(out,1,1,nTr);
                    
                case {'Sine' 'SineWave'}
                    
                    genSin = sin((2*pi*p.Results.Frequency/1e3).*tm + p.Results.PhaseOffset)...
                        .* 0.5 + 0.5;
                    out = (out-p.Results.Baseline) .* genSin + p.Results.Baseline;
                    out = repmat(out,1,1,nTr);
                    
                case 'SmoothPulse'
                    
                    derv = [0 diff(out)];
                    s(1) = p.Results.RiseKernel;
                    s(2) = p.Results.FallKernel;
                    out = cumsum(conv(derv .* (derv>0), normpdf(0:s(1)*8, s(1)*4,s(1)/3), 'same')...
                        + conv(derv .* (derv<0), normpdf(0:s(2)*8, s(2)*4,s(2)/3), 'same'));
                    out = repmat(out,1,1,nTr);
                    
                case {'Noise', 'WhiteNoise'}
                    
                    if p.Results.FixedSeed
                        ns = repmat(rand(1,dur),1,1,nTr);
                    else
                        ns = rand(1,dur,nTr);
                    end
                    
                    out = repmat((out-p.Results.Baseline), 1,1,nTr);
                    out = out .* ns + p.Results.Baseline;
                    

                    
                otherwise
                    
                    error('Type invalid')
            end
            
            
            if obj.isNull
                obj.I = out;
            else
                obj.I(1,:,whichTrial) = out;
            end
            
            if ~isempty(setdiff(1:obj.numTrials, whichTrial))
                assert(all(obj.t == tm), 'Time vector conflict with existing trials');
            else
                obj.t = tm;
            end
            
        end
    end
    methods (Static)
        function out = rngfn(vec, rng)
            % check where array is in range
            out = vec >= rng(1) & vec < rng(2);
        end
    end
end
