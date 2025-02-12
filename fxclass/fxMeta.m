classdef fxMeta
    properties
        Dir char = ''
        ModelSrc char = ''
        ParamSrc char = ''
    end
    methods
        function obj = fxMeta(name, varargin)
            % create object either from pointer (name interpreted as dir in models/),
            % or specifying ModelSrc and ParamSrc directly.
            % Use modelname.m for ModelSrc and getPARAM.m for ParamSrc
            %
            % functions defining model and param need to conserve the
            % following input-output logic:
            %
            % param:    param = ParamSrc()
            %
            % model:    g = ModelSrc(PARAM)
            %
            
            p = inputParser;
            p.addParameter('ModelSrc', []);
            p.addParameter('ParamSrc', []);
            
            p.parse(varargin{:})
            
            if ~isempty(p.Results.ModelSrc) && ~isempty(p.Results.ParamSrc)
                obj.Dir = './';
                obj.ModelSrc = p.Results.ModelSrc;
                obj.ParamSrc = p.Results.ParamSrc;
            else
                obj.Dir = ['models/' name '/'];
                obj.ModelSrc = name;
                obj.ParamSrc = 'getPARAM';
            end

        end
        %%
        function [g, PARAM] = init(obj, varargin)
            obj.mdpath(1)
            PARAM = feval(obj.ParamSrc, varargin{:});
            
            g = feval(obj.ModelSrc, PARAM);
            obj.mdpath(0)
        end
        %%
        function [out, grid] = map(obj, paramrange, varargin)
            % param search using a single iteration for each param
            % combination. Pass up to 4 parameters;
            % paramrange should be of the form {'param1', 0:0.01:1, 'param2', 5:2:17};
            % if more than one param should be varied together, use cell
            % array {{'param1', 'param2'}, [...]}, values will be identical
            % between these parameters.
            % pass args to all proc steps as 
            % {'funcname', {'arg1', val1, arg2, val2, ...}}
            % first func should be 'init' for izMeta.init()
            
            NPass = numel(varargin);
            
            assert(NPass>=1 && strcmpi(varargin{1}{1}, 'init'), 'Pass init args')
            
            
            NParam = numel(paramrange)/2;
            assert(NParam < 5, 'Pass up to 4 parameters');
            
            params = paramrange(1:2:end);
            for iparam = 1:NParam % turn all into cells
               if ~iscell(params{iparam})
                  params{iparam} = {params{iparam}}; 
               end
            end
            
            vals = paramrange(2:2:end);
            vals(NParam+1:4) = {0};
            
            [grd{1}, grd{2}, grd{3}, grd{4}] = ndgrid(vals{:});
            [grdIdx{1}, grdIdx{2}, grdIdx{3}, grdIdx{4}] = ind2sub(size(grd{1}),1:numel(grd{1}));

            if NParam == 1
                grd{1} = vals{1};
            end
            
            for i = 1:4
                grd{i} = grd{i}(:);
            end
            
            N = numel(grd{1}); % num of combinations
            out = {};
%             wb = waitbar(0, 'PARAM MAPPING PROGRESS');
            
            usePar = 1;

            if usePar && isempty(gcp('nocreate'))
                parpool;
            end
            ppm = ParforProgMon(sprintf('Progress: %d Iterations', N), N);
            
            parfor (i = 1:N, usePar)
%                 waitbar(i/N, wb, 'PARAM MAPPING PROGRESS')
                initargs = varargin{1}{2};
                for iParam = 1:NParam
                    for iParamStr = 1:numel(params{iParam})
                        
                        % delete requested parameter from the original
                        % array of args in case it was passed 
                        idx = find(cellfun(@(x) strcmpi(x, params{iParam}{iParamStr}), initargs));
                        if ~isempty(idx)
                           initargs(idx:idx+1) = []; 
                        end
                        % add the value for this iteration
                        initargs{end+1} = params{iParam}{iParamStr};
                        initargs{end+1} = grd{iParam}(i);
                    end
                end
                
                g = obj.init(initargs{:});
                
                tmp = g;
                for k = 2:NPass
                    func = varargin{k}{1};
                    args = varargin{k}{2};
                    
                    if numel(varargin{k}) > 2
                        argPos = varargin{k}{3};
                    else
                        argPos = 1;
                    end
                    
                    nOut = nargout(func);
                    if nOut == -1
                       nOut = nargout(args{1}{1}) 
                    end
                    tmpNOut = cell(1,nOut);
                    [tmpNOut{:}] = feval(func, tmp, args{:});
                    tmp = tmpNOut{argPos};
                end
                out{i} = tmp;  
                ppm.increment()
            end
            
            delete(ppm);
            
%             out = [out{:}];
            grid = [];
            for iParam = 1:NParam
                for iParamStr = 1:numel(params{iParam})
                    me = params{iParam}{iParamStr};
                    grid.(me) = grd{iParam};
                end
                grid.idx(iParam,:) = grdIdx{iParam};
            end
            
        end
        %%
        function PARAM = getDefaultParam(obj)
            obj.mdpath(1)
            PARAM = feval(obj.ParamSrc);
            obj.mdpath(0)
        end
        %%
        function mdpath(obj, flag)
            % put here in case we want to expand pointers
            if flag
                addpath(obj.Dir);
            else
                rmpath(obj.Dir);
            end
        end
    end
end