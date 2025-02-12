function param = pvn_getNodeHandlesDualStream(param)
% converts the node handles and labels created by pvn_param to the
% dual-stream model in which every node exists as _L and _R. New handles
% are unique to each hemisphere! (but still shared between areas of the
% same hemisphere).
% Handle to IGX spiking neuron is currently still shared between
% hemispheres.

hemi = {'L' 'R'};


lbl = cellfun(@(h) cellfun(@(x) [x '_' h], param.NodeLabels, 'UniformOutput', false), hemi, 'UniformOutput', false);
param.NodeLabels = [lbl{:}];

param.Nodes = [param.Nodes param.Nodes.copy];
param.nodeInArea = repmat(param.nodeInArea,1,2);

% Syn handles are created by pvn_connect which handles both single-stream
% and dual-stream 

end