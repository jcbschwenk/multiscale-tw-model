function [g, param, stim, prior] = pvn_mfmodel_dual(param)

% double hemispheric laminar model with internal alpha generators in L5/6

if ~nargin
    param = pvn_param();
end

if nargin < 2
    name = ['PVN-' genRandCode(8)];
end

%% ARCHITECTURE:
g = fxGraph();

g.FI = param.fiCurve;
g.spkIConv = param.spkIConv;

isInArea = @(nd,i) any(strcmpi(param.NodeLabels, nd)) && ismember(i, param.nodeInArea{strcmpi(param.NodeLabels, nd)});

for iArea = param.Areas
    for iType = param.NodeLabels
        ndType = char(iType);
        ndHandle = param.Nodes(strcmpi(param.NodeLabels, ndType));
        if isInArea(ndType, iArea)
            g.addNode(ndHandle, [ndType num2str(iArea)]);
        end
    end
end


%% ADD EXTERNAL CURRENTS:

[g, param, stim, prior] = pvn_addExtCurr(g, param);


%% CONNECTIONS:

[g, param] = pvn_connect(g, param);

%%
g.name = name;

end
%%
function out = genRandCode(len, seed)

if nargin > 1
    rnd(seed);
end

symbols = ['a':'z' '0':'9'];
nums = randi(numel(symbols),[1 len]);
out = symbols (nums);
end