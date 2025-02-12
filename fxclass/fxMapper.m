classdef fxMapper < matlab.mixin.Copyable
    properties
        Graph(1,1) fxGraph
        ObjHandles(1,:) cell % all below are 1 x NParams cell arrays
        Param(1,:) cell % 
        Vals(1,:) cell %
        Dim(1,:) double % 1 x NParams
        StackFn(1,:) cell = {@(x, args) run(x, args{:})}
        StackLabel(1,:) cell = {'run'}
        StackVarArgIn(1,:) cell = {}
        StackPosArgIn(2,:) double = [0; 1] % first row is fn, second is out pos
        StackPosArgOut(1,:) double = 1 % which output arg to return
    end
    properties (Access = private)
        InitVals(1,:) double
    end
    %%
    methods
        function obj = fxMapper(varargin)
            isGraph = cellfun(@(x) isa(x, 'fxGraph'), varargin);
            if ~isempty(isGraph)
                obj.Graph = varargin{isGraph};
            end
            
            other = find(cellfun(@(x) iscell(x), varargin));
            for i = 1:numel(other)
                args = varargin(other(i));
                obj.add(args{:});
            end
        end
        %%
        function obj = add(obj, varargin)
            % add grid variables using the following structure:
            % {Object, ParamName, ParamVals [, Dimension]}.
            % Object has to be present in obj.Graph and ParamName has to be
            % a parameter of Object. Size of ParamVals array has to match
            % other definition already added to Dimension.
            def = varargin;
            if numel(def) < 4
                def{4} = obj.NumDims + 1;
            end
            [targObj, targParam, paramVals, targDim] = deal(def{:});
            
            assert(obj.checkGraphContains(targObj), 'Target object not in graph');
            assert(isprop(targObj, targParam), 'Target parameter not a property of the object');
            
            if obj.NumDims > 0 && targDim <= obj.NumDims
                [~, dimSz] = obj.checkDim(targDim);
                assert(dimSz == numel(paramVals), 'Dimension size mismatch with existing parameter');
            end
   
            obj.ObjHandles{end+1} = targObj;
            obj.Param{end+1} = targParam;
            obj.Vals{end+1} = paramVals;
            obj.Dim(end+1) = targDim;
            
        end
        %%
        function obj = clearGrid(obj)
            obj.ObjHandles = {};
            obj.Param = {};
            obj.Vals = {};
            obj.Dim = [];
        end
        %%
        function [obj, k] = addToStack(obj, varargin)
            % add new function to the stack. Input arguments are positional
            % and must be: 1) label 2) function handle. 3) additional arguments (any
            % after first) to the function as cell array 4) input/output
            % arg position: [2x1 double], first element is index of function
            % from which to take the output as input, second element is the
            % position. Default is [k-1; 1] for function k (i.e. using
            % first output of previous function in the stack). 5) output
            % arg pos to return (as output of mapper.run()), default is 1.
            
            obj.StackFn{end+1} = varargin{2};
            k = numel(obj.StackFn);
            
            obj.StackLabel{k} = varargin{1};
            

                  
            obj.StackVarArgIn{k} = {};
            if nargin > 3
                obj.StackVarArgIn{k} = varargin{3};
            end
            
            obj.StackPosArgIn(:,k) = [k-1; 1];
            if nargin > 4
                obj.StackPosArgIn(1:numel(varargin{4}),k) = varargin{4};
            end
            
            obj.StackPosArgOut(k) = 1;
            if nargin > 5
                obj.StackPosArgOut(k) = varargin{5};
            end
            
        end
        %%
        function obj = clearStack(obj)
            obj.StackFn(2:end) = [];
            obj.StackLabel = {'run'};
            obj.StackVarArgIn = {};
            obj.StackPosArgIn = [0; 1];
            obj.StackPosArgOut = 1;
        end
        %%
        function [varargout] = run(obj, varargin)
            % run the function stack on the specified parameter grid.
            % Output arguments are in the following order:
            % 1:  grid definition (contains param names and values)
            % 2 to N: output of function stack in reverse order (i.e.
            % varargout{end} will be output of graphObj.run).
            % Alternatively, if 'struct' is passed, return arg is struct
            % with named fields containing all output.
            
            if ismember('struct', varargin)
                returnStruct = true;
                varargin(strcmpi('struct', varargin)) = [];
            else
                returnStruct = false;
            end
            
            [grd, idxGrd] = obj.getGrid();
            
            N = size(grd,2); % this is the number of runs
            
            
            K = numel(obj.StackFn); 
            
            varargout = cell(1,K);
            
            if nargin > 1 && ~isempty(varargin)
                obj.StackVarArgIn(1) = varargin;
            else
                obj.StackVarArgIn(1) = {{}};
            end
            
            tmpVarOut = cell(1,K);
            for k = 1:K
                usedBy = find(obj.StackPosArgIn(1,:) == k);
                if isempty(usedBy)
                    nOut = obj.StackPosArgOut(k);
                else
                    nOut = max(obj.StackPosArgOut(k), max(obj.StackPosArgIn(2, usedBy)));
                end
                tmpVarOut{k} = cell(1, nOut);
            end
            
            obj.getInitVals(); % fetch initial values of all objects
            
            wb = waitbar(0, sprintf('Grid run [%d D] [%d iterations]', obj.NumDims, N));
            for iRun = 1:N
                
                waitbar(iRun/N, wb, sprintf('Grid run [%d D] [%d iterations], iteration %d/%d', obj.NumDims, N, iRun, N));
                
                obj.setVals(idxGrd(:,iRun)) % set to values for this run
                
                posArg = obj.Graph; % first input arg is always the graph
                
                for kFun = 1:K
                    if kFun > 1
                        if obj.StackPosArgIn(1,kFun) == 0
                            posArg = obj.Graph;
                        else
                            posArg = tmpVarOut{obj.StackPosArgIn(1,kFun)}{obj.StackPosArgIn(2,kFun)};
                        end
                    end
                    varArgs = obj.StackVarArgIn(kFun);
                    [tmpVarOut{kFun}{:}] = obj.StackFn{kFun}(posArg, varArgs{:}); % function call
                    varargout{K-kFun+1}{end+1} = tmpVarOut{kFun}{obj.StackPosArgOut(kFun)}; % collect output
                end
            end
            delete(wb)
            obj.resetToInitVals(); % reset all to initial values

            % varargout is now 1 x K cell array in reverse order of
            % function call stack
            

            grdDef = [];
            grdDef.val = grd;
            grdDef.idx = idxGrd;
            grdDef.param = obj.Param;
            
            
            
            if returnStruct
                out = [];
                out.grd = grdDef;
                varargout = varargout(end:-1:1);
                for iFn = 1:numel(varargout)
                   out.(obj.StackLabel{iFn}) = varargout{iFn};
                   try
                       out.(obj.StackLabel{iFn}) = [out.(obj.StackLabel{iFn}){:}];
                   catch me
                       warning('Leaving output of function #%d as cell array:', iFn)
                       if strcmpi(me.identifier, 'MATLAB:catenate:structFieldBad')
                           warning('Struct fields for different runs do not match, likely due to a default empty output')
                       end
                   end
                end
                varargout = {out};
            else
                varargout = [{grdDef} varargout];
            end
            
        end
        %%
        function [grd, idxGrd] = getGrid(obj)
            for i = 1:obj.NumDims
               [~, sz(i)] = obj.checkDim(i);
            end
            idxGridList = arrayfun(@(x) 1:sz(x),1:obj.NumDims, 'UniformOutput', false);
            
            idxGrd = cell(1,obj.NumDims);
            [idxGrd{:}] = ndgrid(idxGridList{:});
            
            idxGrd = cellfun(@(x) x(:)', idxGrd, 'UniformOutput', false); % flatten
            
            N = numel(obj.Param); 
            grd = arrayfun(@(x) obj.Vals{x}(idxGrd{obj.Dim(x)}),...
                1:N, 'UniformOutput', false); % index into the actual values
            
            grd = cat(1, grd{:});
            idxGrd = cat(1, idxGrd{:});
            
        end
        %%
        function [out, sz] = checkDim(obj, whichDim)
            
            if nargin < 2
                for i = unique(obj.Dim)
                    [out(i), sz(i)] = obj.checkDim(i);
                end
                return
            end
            n = cellfun(@(x) numel(x), obj.Vals(obj.Dim == whichDim));
            out = numel(unique(n)) == 1;
            sz = out * n(1);
        end
        %%
        function out = checkGraphContains(obj, targObjHandle)
            
            inNodes = any(obj.Graph.Nodes == targObjHandle);
            inCurr = any(obj.Graph.Curr == targObjHandle);
            inSyn = any(obj.Graph.Syn(:) == targObjHandle);
            allSpk = [obj.Graph.Nodes(obj.Graph.Nodes.hasSpk).SpkChild];
            inSpk = any(allSpk(:) == targObjHandle);
            
            out = inNodes || inCurr || inSyn || inSpk;
            
        end
        %%
        function n = NumDims(obj)
            n = max(obj.Dim);
            if isempty(n)
               n = 0; 
            end
        end
        %%
    end
    %%
    methods (Access = private)
        %%
        function obj = getInitVals(obj)
            for i = 1:numel(obj.Param)
                obj.InitVals(i) = obj.ObjHandles{i}.(obj.Param{i});
            end
        end
        %%
        function obj = resetToInitVals(obj)
            for i = 1:numel(obj.Param)
                obj.ObjHandles{i}.(obj.Param{i}) = obj.InitVals(i);
            end
        end
        %%
        function obj = setVals(obj, idx)
            % set param vals to the values at pos idx
            for i = 1:numel(obj.Param)
                obj.ObjHandles{i}.(obj.Param{i}) = obj.Vals{i}(idx(obj.Dim(i)));
            end
        end
    end
    %%
    methods (Static)
        %%
        function [isType, typeStr] = checkType(in)
            validTypes = {'fxNode' 'fxSpkNeurIzh' 'fxSyn' 'fxCurr'};
            isType = find(cellfun(@(x) isa(in, x), validTypes));
            if isempty(isType)
               isType = 0;
               typeStr = '';
               return
            end
            typeStr = validTypes{isType};                
        end
    end
end