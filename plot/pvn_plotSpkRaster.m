function pvn_plotSpkRaster(in, g, varargin)

p = inputParser;
p.addParameter('Node', {});
p.addParameter('Dark', false)

parse(p, varargin{:});

nodes = p.Results.Node;

if isempty(nodes)
    nodes = in.spkNodeIdx;
else
    nodes = g.getNodeIdx(nodes);
end
nodes = nodes(nodes > 0);

if p.Results.Dark
    fgcol = [1 1 1]; bgcol = [0 0 0];
else
    fgcol = [0 0 0]; bgcol = [1 1 1];
end

neur2NodeSel = in.neur2Node(ismember(in.neur2Node, nodes));

dur = numel(in.t);
nNeur = numel(neur2NodeSel);

[~, idx] = ismember(nodes, in.spkNodeIdx);

assert(any(idx), 'Requested nodes are not spiking nodes')

nodes(idx == 0) = [];
idx(idx == 0) = [];

spkSel = in.spk(idx);
allSpk = cat(1, spkSel{:});

mmt = zeros(nNeur,dur);
if numel(nodes) == 1 % single node
    yVals = allSpk(:,3)'; % index in node
elseif numel(nodes) == numel(in.spkNodeIdx) % all nodes
    yVals = allSpk(:,2)'; % index in full graph
else % node selection
    yVals = allSpk(:,3)';
    k = 0;
    for iNode = 1:numel(nodes)
        thisIdx = in.neur2Node(allSpk(:,2)) == nodes(iNode);
        yVals(thisIdx) = yVals(thisIdx) + k;
        k = k + max(yVals(thisIdx));
    end
end

xVals = allSpk(:,1)';

mmt(sub2ind(size(mmt), yVals, xVals)) = 1;

raster = imagesc(1:dur, 1:nNeur, mmt);
raster.AlphaData = mmt;
set(gca, 'Colormap', fgcol);
set(gca, 'Color', bgcol)
caxis([0 1])

txtYPos = [1 cumsum(cellfun(@(x) numel(unique(x(:,2))), in.spk(idx(2:end))))];
txtYPos = txtYPos + 5; % offset
for i = 1:numel(idx)
    text(min(xlim), txtYPos(i), g.Labels{nodes(i)},...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top',...
        'FontWeight', 'bold');
end

end