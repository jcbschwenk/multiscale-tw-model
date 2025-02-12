classdef fxSyn < matlab.mixin.Copyable
    properties
        W double = 0 % weight
        T double = 0 % delay
        Tau double = 1 % synaptic decay parameter
        G double = 0 % saturation parameter gamma
    end
    methods
        function obj = fxSyn(varargin)
            
            p = inputParser;
            for iProp = cell(properties('fxSyn'))'
                p.addParameter(char(iProp), []);
            end
            
            p.parse(varargin{:})
            
            for iProp = cell(properties(obj))'
                prop = char(iProp);
                if ~isempty(p.Results.(prop))
                    obj.(prop) = p.Results.(prop);
                end
            end
        end
        %%
        function [W, T, Tau, G] = getParam(obj)
            [W, T, Tau, G] = deal(zeros(size(obj)));
            W(:) = [obj.W];
            T(:) = [obj.T];
            Tau(:) = [obj.Tau];
            G(:) = [obj.G];
        end
    end
end
