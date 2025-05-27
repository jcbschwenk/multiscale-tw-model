% Parameter search for feed-up weight w/ eeg projection, in a
% single-stream model.

function [f, f2] = pvn_fdupWeightSingleStream()

basefolder = pwd; % set manually
tmpfolder = [basefolder 'fdupWeight1D/']; mkdir(tmpfolder)
%%

w = 0:0.1:2.2;

% sum same delay:
pulThru(2) = fUpWeight_noFBkFUp(w, 0);


%% PLOTS


delays = {'Double' 'Same'};

iDel = 2;

f = figure;
tiledlayout(3,1);
pfw(:,:,iDel) = cat(2,pulThru(iDel).dst.pfw);
pbw(:,:,iDel) = cat(2,pulThru(iDel).dst.pbw);

curr = pulThru(iDel).dst.curr;
% -> now a matrix of nStepsCurr x nStepsWeight x delay

nexttile([2 1])
contourf(w, curr, pfw(:,:,iDel), 'EdgeColor', 'none');
hold on
plot(w, [pulThru(iDel).thresh], 'k:', 'LineWidth', 2)
caxis([0 1])
colormap(brewermap(64, 'RdBu'))
cb = colorbar;
cb.Ticks = [0 0.5 1];
title(cb, 'P(FW)')
xlim([min(w) max(w)])

xticks([])
ylabel('I_{STIM} [nA]')

nexttile([1 1])
plot(w, log(pfw(1,:,iDel)./pbw(1,:,iDel)), 'k-', 'LineWidth', 2);
xlabel('Weight Pul --> Cx')
text(min(xlim)+0.05*diff(xlim), 0.95* max(ylim), 'log(P(FW | I_{STIM} = 0)/P(BW | I_{Stim} = 0))',...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold')
ylim([min(ylim) 0.5]);
xlim([min(w) max(w)])
yline(0, 'k-')
sgtitle([delays{iDel} ' Delay'])


set(gcf, 'Color', 'w')

%%
f2 = figure;
dirlabel = {'pfw' 'pbw'};
for idir = 1:2
subplot(1,2,idir)

dstMin = pulThru(iDel).dst(1);
[~, idx] = min(abs(w-1.5));
dstMax = pulThru(iDel).dst(idx);

hold on
ln(1) = pvn_shplot(dstMin.curr, dstMin.(dirlabel{idir}), dstMin.([dirlabel{idir} 'CI']), 'LineWidth', 2);
ln(2) = pvn_shplot(dstMax.curr, dstMax.(dirlabel{idir}), dstMax.([dirlabel{idir} 'CI']), 'LineWidth', 2);
xlabel('Input Current [nA]')
ylabel('P(FW)')

xlim([min(dstMax.curr) max(dstMax.curr)] + [-0.1 0.1])
xline(min(dstMax.curr), 'k');
xline(max(dstMax.curr), 'k');
ylim([0 1.1])
legend(ln, {'Cx' 'Cx+Pul (w=1)'}, 'Location', 'NorthWest')
legend('boxoff')
end
set(gcf, 'Color', 'w')


%%
save([tmpfolder 'FdUpWeight.mat'])
%%
end
%%
function [g, param, m] = init()

[g, param] = pvn_mfmodel();
m = fxMapper(g);

end
%%
function getStack(m, param)
% Set stack as full dynamic state threshold estimation (runs on graph, i.e.
% replaces 'run' as first fn in the stack):

dstArgs = {...
    param,...
    'MaxCurr', 1.5,...
    'PriorScale', 0.3,...
    'NumSteps', 15,... 
    'NumTrialsPerStep', 10,...
    'Figure', false,...
    'SaveOutput', false};

m.StackFn{1} = @(x, args) pvn_dynStateThresh(x, dstArgs{:}, args{:});
m.StackLabel{1} = 'dst';

% extract threshold and stim-off-baseline (purely for convenience):
m.addToStack('thresh',...
    @(x, args)...
    x.thresh,...
    {}, [1; 1]);
m.addToStack('lower',...
    @(x, args)...
    x.lower,...
        {}, [1; 1]);
end
%%
function out = fUpWeight_noFBkFUp(w, isDoubleDelay)

[g, param, m] = init();

getStack(m, param);

fUpSyn = param.Syn(strcmpi(param.SynLabels, 'FUp_Pul_L4X'));
fDnSyn = param.Syn(strcmpi(param.SynLabels, 'FDn_IGX_Pul'));

m.add(fUpSyn, 'W', w);

% fix delays at same or double (independent of default):

if isDoubleDelay
    m.add(fUpSyn, 'T', param.conn.interAreaDelay, 2);
    m.add(fDnSyn, 'T', param.conn.interAreaDelay, 2);
else
    m.add(fUpSyn, 'T', ceil(param.conn.interAreaDelay/2), 2);
    m.add(fDnSyn, 'T', floor(param.conn.interAreaDelay/2), 2);
end

out = m.run('struct');

end
