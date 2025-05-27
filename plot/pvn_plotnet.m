function pvn_plotnet(g, param)
% plots the connectivity structure of the final graph 

isdual = any(contains(param.NodeLabels, '_R'));

%%
specs.ILSt = {[1 3], 'k'};
specs.ILStIN = {[3 2], 'r'};

specs.ILPr = {[1 5], 'k'};
specs.ILPrIN = {[3 6], 'r'};

specs.L4X = {[1 4], 'k'};
specs.L4IN = {[1.3 4.3], 'r'};
specs.L4FF = {[3 4], 'g'};

specs.SGX = {[2 2], 'r'};
specs.SGE = {[3 2], 'y'};
specs.SGIN = {[3 3], 'g'};
specs.IGX = {[2 5], 'b'};
specs.IGIN = {[2 6], 'c'};
specs.Pul = {[3 9], 'k'};
specs.PulIN = {[4 9], 'r'};

%%

cxwidth = 3;
cxoffset = 2;

%%

% figure
xlim([-6 (param.N+1)*(cxwidth+cxoffset)])
ylim([-1 11]);
set(gca, 'ydir', 'reverse');
hold on


%%
if ~isdual
    nds = g.Labels;
    tps = param.NodeLabels;
else
    whichHemi = 'L';
    nds = g.Labels(contains(g.Labels, ['_' whichHemi]));
    tps = param.NodeLabels(contains(param.NodeLabels, ['_' whichHemi]));
    
    specs = cell2struct(struct2cell(specs), cellfun(@(x) [x '_' whichHemi], fields(specs), 'UniformOutput', false));
end


NNds = numel(tps);
[idx, xy] = deal([]);
[lbl, tp] = deal({});

for i = 0:param.N+1
    for iType = 1:NNds
        mtp = tps{iType};
        
        if ~ismember(i, param.nodeInArea{strcmpi(param.NodeLabels, mtp)})
            continue
        end
        idx(end+1) = g.getNodeIdx([mtp num2str(i)]);
        tp{end+1} = mtp;
        lbl{end+1} = [mtp num2str(i)];
        xy(end+1,:) = specs.(mtp){1} + [(i-1) * (cxwidth+cxoffset) 0];
        scatter(xy(end,1), xy(end,2), [], specs.(mtp){2}, 'filled');         
    end
    
    cxXLim = (i-1)*(cxwidth+cxoffset) + [0 cxwidth+1];
    
    if i == 0 || i == param.N+1
        cxYLim = [1 3] + (i==param.N+1)*3 + 0.5;

        line([cxXLim(1) .* [1 1] cxXLim(2) .* [1 1] cxXLim(1)],...
            [cxYLim(2) cxYLim(1) .* [1 1] cxYLim(2) .* [1 1]],...
            'LineStyle', '-', 'Color', 'k')
      
        text(mean(cxXLim), min(cxYLim)-0.5, 'IL', 'Color', 'k', 'FontWeight', 'bold',...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    else
            cxYLim = [1 6] + 0.5;

        line([cxXLim(1) .* [1 1] cxXLim(2) .* [1 1] cxXLim(1)],...
            [cxYLim(2) cxYLim(1) .* [1 1] cxYLim(2) .* [1 1]],...
            'LineStyle', '-', 'Color', 'k')
        line(cxXLim, [4 4]-0.5,'LineStyle', ':', 'Color', 'k')
        line(cxXLim, [5 5]-0.5, 'LineStyle', ':', 'Color', 'k')
        
        text(mean(cxXLim), min(cxYLim)-0.5, ['Cx' num2str(i)], 'Color', 'k', 'FontWeight', 'bold',...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    
end




%%

for ind = 1:numel(idx)
    [~,~,rec] = g.getRecipients(lbl{ind});
    
    for irec = 1:numel(rec)
        other = find(idx == rec(irec));
        lptsX = [xy(ind,1) xy(other,1)];
        lptsY = [xy(ind,2) xy(other,2)];
        
        rnd = rand/2+0.25;
        
        ww = g.W;
        ww = ww(idx(ind),rec(irec));
        
        text(sum([rnd 1-rnd].*lptsX), sum([rnd 1-rnd].*lptsY),...
            sprintf('%1.1f', ww), 'Color', specs.(tp{ind}){2});
        
        if ww > 0 
            prop = 0.9;
            scatter(sum([1-prop prop].*lptsX), sum([1-prop prop].*lptsY), 40, specs.(tp{ind}){2}, 'filled')
                line(lptsX(1)+[0 prop].*diff(lptsX),...
                    lptsY(1)+[0 prop].*diff(lptsY),...
                    'Color', specs.(tp{ind}){2}, 'LineWidth', 1.5);

        elseif ww < 0 
              prop = 0.75;
            scatter(sum([1-prop prop].*lptsX), sum([1-prop prop].*lptsY), 80, specs.(tp{ind}){2}, '_',...
                'LineWidth', 2)
            line(lptsX(1)+[0 prop].*diff(lptsX),...
                    lptsY(1)+[0 prop].*diff(lptsY),...
                    'Color', specs.(tp{ind}){2}, 'LineWidth', 1.5);

        end
    end
end
%%
axis off
end