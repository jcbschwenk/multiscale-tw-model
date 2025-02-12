classdef fxSpkNeurIzh < matlab.mixin.Copyable
    properties
        % Izh parameters a,b,c,d
        a double = 0.02 
        b double = 0.2  
        c double = -65  
        d double  = 8 
        e(2,1) double = [2 5] % gaussian noise [mu sigma]
    end
    methods
        function obj = fxSpkNeurIzh(varargin)
            
            p = inputParser;
            for iProp = cell(properties('fxSpkNeurIzh'))'
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
        
    end
end
