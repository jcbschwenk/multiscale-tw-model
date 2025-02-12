classdef fxNode < matlab.mixin.Copyable
    properties
        lm double = 270 % FI param lambda (excitability), in Hz/nA
        b double = 108 % FI parameter b, in Hz
        c double = 0.154 % FI parameter c, in sec 
        Ib double  = 0.3 % base current, in nA
        e double = 0.004 % noise sigma
        eTD double = 3 % decay tau for noise in ms
        SpkChild fxSpkNeurIzh 
    end
    methods
        function obj = fxNode(varargin)
            
            p = inputParser;
            
            props = cell(properties('fxNode'))';
            
            for iProp = props
                p.addParameter(char(iProp), []);
            end
            
            p.parse(varargin{:})
            
            for iProp = props
                prop = char(iProp);
                if ~isempty(p.Results.(prop))
                    obj.(prop) = p.Results.(prop);
                end
            end
        end
        %%
        function obj = addSpkChildren(obj, N, varargin)
            % adds N neurons as spk children. All parameters in fxSpkNeurIzh
            % must be passed as fixed or a range of [min max]. Pass args as
            % name value pairs.
            
            p = inputParser();
            
            props = cell(properties('fxSpkNeurIzh'))';
            
            for iProp = props
                p.addParameter(char(iProp), []);
            end
            
            p.parse(varargin{:})
            obj.SpkChild = fxSpkNeurIzh.empty(0,N);
            [obj.SpkChild.e] = deal(p.Results.e);
            
            for iProp = setdiff(props, 'e')
                val = p.Results.(char(iProp));
                
                val = [min(val) max(val)];
                val = rand(1,N).*diff(val) + val(1);
                for iNeur = 1:N
                    obj.SpkChild(iNeur).(char(iProp)) = val(iNeur);
                end
            end
            
        end
        %%
        function out = hasSpk(obj)
           % returns true if node has any non-empty spiking children
           out = false(1,numel(obj));
           for i = 1:numel(obj)
              out(i) = ~isempty(obj(i).SpkChild);
           end
        end
        %%
        function out = numSpkChildren(obj)
            out = zeros(1,numel(obj));
            for i = 1:numel(obj)
                out(i) = numel(obj(i).SpkChild);
            end
        end
        %%
        function [a,b,c,d,e,ID] = getSpkChildParam(obj)
           % for a 1D array of nodes, get izh parameters for all spk
           % children with node ID as a separate output arg
                    
           assert(ismatrix(obj) && any(size(obj) == 1), 'Input must be single-row or -column node vector.');
           
           if ~any(obj.hasSpk)
              [a,b,c,d,e,ID] = deal([]);
              return
           end
           
           allSpk = [obj(obj.hasSpk).SpkChild];
           a = [allSpk.a];
           b = [allSpk.b];
           c = [allSpk.c];
           d = [allSpk.d];
           e = [allSpk.e];
           
           numSpk = obj.numSpkChildren;
           ID = arrayfun(@(x,y) repmat(x,1,y), find(numSpk > 0), numSpk(numSpk > 0), 'UniformOutput', false);
           ID = [ID{:}];
        end
        %% 
    end
end
