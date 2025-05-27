function pvn_plotMap(in, g, param, ndTypes)

if nargin < 4
    ndTypes = param.nodeTypes;
end

if ~iscell(ndTypes)
   ndTypes = {ndTypes}; 
end

k = 1;
for iType = 1:numel(ndTypes)
    nd = ndTypes{iType};
    
    for i = 0:param.N+1
        ndIdx = g.getNodeIdx([nd num2str(i)]);
        if ~ndIdx
            continue
        end
        map(k,:) = in.r(ndIdx,:);
        ndlbl(k) = iType;
        k = k+1;
    end
    
    
end

imagesc(in.t, 1:size(map,1), map)

sepLineStyle = 'w:';

for iType = 1:numel(ndTypes)
    nd = ndTypes{iType};
    idx = find(ndlbl == iType, 1, 'first');
    if isempty(idx)
        continue
    end
    yline(idx - 0.5, sepLineStyle)
    text(min(xlim)+20, idx-0.5, nd, 'FontWeight', 'bold', 'Color', 'w', 'FontSize', 12,...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end


if any(xlim < 0)
    xline(0, 'w:', 'LineWidth', 2)
end
yticks([])
xlabel('Time [ms]')

if any(caxis > 500)
   caxis([min(caxis) 500]); 
end

set(gcf, 'Color', 'w')
set(gca, 'FontName', 'Arial')

%%
sgtitle(g.name)
end