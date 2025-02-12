function g = pvn_nodePairLayerConnect(g, param, connDir, fromNode, toNode, varargin)
% connects all nodes matching fromNode and toNode based on numbering,
% following the rule defined by connDir. connDir is an array of target
% layers as integers relative to sending layer 0 (e.g. [-1 1] for symmetric
% feedforward and feedback). Pass connection args after node labels


fromNodeIdx = strcmpi(param.NodeLabels, fromNode);
toNodeIdx = strcmpi(param.NodeLabels, toNode);

assert(any(fromNodeIdx) && any(toNodeIdx), 'Bad node reference')

fromAreas = param.nodeInArea{fromNodeIdx};
toAreas = param.nodeInArea{toNodeIdx};

areaDist = toAreas - fromAreas';

for targDist = unique(connDir)
    [fromIdx, toIdx] = find(areaDist == targDist) ;
    fromIdx = fromAreas(fromIdx);
    toIdx = toAreas(toIdx);
    
    for k = 1:numel(fromIdx)
        
        fromNodeLabel = [fromNode num2str(fromIdx(k))];
        toNodeLabel = [toNode num2str(toIdx(k))];
        
        fprintf('\nConnecting %s to %s\n', fromNodeLabel, toNodeLabel);
        
        g.connect(fromNodeLabel, toNodeLabel, varargin{:});

    end
end

end