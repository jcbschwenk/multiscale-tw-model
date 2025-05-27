% Parameter search for feed-up weight w/ eeg projection, in a
% dual-stream model where the two hemispheres are not inter-connected.

function [f1, f2] = pvn_fdupWeightDualStreamParallel()

basefolder = pwd; % set manually
tmpfolder = [basefolder 'fdupWeight2DParallel/']; mkdir(tmpfolder)
%%
w = 0:0.1:2.2;

fixedWeightRightHemi = 0.2;

pulThru(1) = fUpWeight(w, fixedWeightRightHemi, 1);

%% PLOTS
dels = {'Same' 'Double'};

hemis = {'L' 'R'};


iDel = 1;

%%
f1 = figure;

tiledlayout(1,5);

for iHemi = 1:2
    
    
    tmp = arrayfun(@(x) x.pfw(:,iHemi), pulThru(iDel).dst, 'UniformOutput', false);
    pfw(:,:,iHemi) = [tmp{:}];
    curr = pulThru(iDel).dst(1).curr;
    % -> now a matrix of nStepsCurr x nStepsWeight x delay
    
    nexttile((iHemi-1)*2 + 1)
    contourf(w, curr, pfw(:,:,iHemi), 'EdgeColor', 'none');
    hold on
    plot(w, [pulThru(iDel).(['thresh_' hemis{iHemi}])], 'k:', 'LineWidth', 2)
    caxis([0 1])
    colormap(brewermap(64, 'RdBu'))
    cb = colorbar;
    cb.Ticks = [0 0.5 1];
    title(cb, 'P(FW)')
    
    xticks([])
    ylabel('I_{STIM} [nA]')
    
    nexttile((iHemi-1)*2 + 2)
    plot(w, [pulThru(iDel).(['lower_' hemis{iHemi}])], 'LineWidth', 2, 'Color', pvn_figCols(['fw' hemis{iHemi}]));
    hold on
    plot(w, [pulThru(iDel).(['upper_' hemis{iHemi}])], 'k-', 'LineWidth', 2, 'Color', pvn_figCols(['fw' hemis{iHemi}]));
    xlabel('Weight Pul --> Cx')
    text(min(xlim)+0.05*diff(xlim), 0.95* max(ylim), 'P(FW | I_{STIM} = 0)',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold')
    
    sgtitle([dels{iDel} ' Delay'])
    
end

nexttile(5)
plot(w, pulThru(iDel).lower_bias, 'm', 'LineWidth', 2);
hold on
plot(w, pulThru(iDel).upper_bias, 'g', 'LineWidth', 2);
set(gcf, 'Color', 'w')

%%
f2 = figure;
tiledlayout(2,3);
[~, idx1] = min(abs(w-0.2));
dstMin = pulThru(iDel).dst(idx1);
[~, idx2] = min(abs(w-2));
dstMax = pulThru(iDel).dst(idx2);

hemis = {'L' 'R'};
dirs = {'fw' 'bw'};

for iDir = 1:2
    for iROI = 1:2
        nexttile((iDir-1)*2+1); hold on
        pvn_shplot(dstMin.curr, dstMin.(['p' dirs{iDir}])(:,iROI), dstMin.(['p' dirs{iDir} 'CI'])(:,:,iROI), 'LineWidth', 2, 'Color', pvn_figCols(['fw' hemis{iROI}]));
        nexttile((iDir-1)*2+2); hold on
        [~, ln(iROI + (iDir-1)*2)] = pvn_shplot(dstMax.curr, dstMax.(['p' dirs{iDir}])(:,iROI), dstMax.(['p' dirs{iDir} 'CI'])(:,:,iROI), 'LineWidth', 2, 'Color', pvn_figCols(['fw' hemis{iROI}]));
    end
    
    for iLim = 1:2
        nexttile((iDir-1)*2+iLim)
        xlim([min(dstMax.curr) max(dstMax.curr)] + [-0.1 0.1])
        xline(min(dstMax.curr), 'k');   
        xline(max(dstMax.curr), 'k');
        ylim([0 1])
    end
end

legend(ln, repmat(hemis, 1, 2), 'Location', 'NorthWestOutside', 'AutoUpdate', 'off')
ax = gca;
ax.YAxis.Visible = 'off';
title(sprintf('W_{Pul->Cx} = %d', w(idx2)))

nexttile(1)
xlabel('Input Current [nA]')
ylabel('P(FW)')
yticks(0:0.2:1)
title(sprintf('W_{Pul->Cx} = %d', w(idx1)))


set(gcf, 'Color', 'w')

%%
save([tmpfolder 'FdUpWeight.mat'])
savefig(f1, [tmpfolder 'FdUpWeight_FIG'])
%%

end
%%
function [g, param, m] = init()

[g, param] = pvn_mfmodel('dual');
m = fxMapper(g);

end
%%
function getStack(m, param)
% Define pipeline for parameter searches:

dstArgs = {...
    param,...
    'MaxCurr', 1.5,...
    'PriorScale', 0.3,...
    'CurrLabel', 'Stim_M',...
    'NumSteps', 15,... % 15
    'NumTrialsPerStep', 20,...
    'Figure', false,...
    'SaveOutput', false};

m.StackFn{1} = @(x, args) pvn_dynStateThresh(x, dstArgs{:}, args{:});
m.StackLabel{1} = 'dst';

hemis = {'L' 'R'};

for iHemi = 1:2
    % extract threshold and stim-off-baseline (purely for convenience):
    m.addToStack(['thresh_' hemis{iHemi}],...
        @(x, args)...
        x.thresh(iHemi),...
        {}, [1; 1]);
    m.addToStack(['lower_' hemis{iHemi}],...
        @(x, args)...
        x.lower(iHemi),...
        {}, [1; 1]);
    m.addToStack(['upper_' hemis{iHemi}],...
        @(x, args)...
        x.upper(iHemi),...
        {}, [1; 1]);
end

end
%%
function out = fUpWeight(w, fixedWeightRightHemi, fixSameDelay)

[g, param, m] = init();

getStack(m, param);

% feed-up from Pulvinar to Cortex, varying left hemisphere only:
idxL = strcmpi(param.SynLabels, 'FUp_Pul_L_L4X_L');
synL = param.Syn(idxL); % this is the target synapse
m.add(synL, 'W', w);

% right hemisphere is fixed:
idxR = strcmpi(param.SynLabels, 'FUp_Pul_R_L4X_R');
synR = param.Syn(idxR);
m.add(synR, 'W', fixedWeightRightHemi);

if fixSameDelay
    
    fDnL = param.Syn(strcmpi(param.SynLabels, 'FDn_IGX_L_Pul_L'));
    fDnR = param.Syn(strcmpi(param.SynLabels, 'FDn_IGX_R_Pul_R'));
    
    dnDelay = floor(param.conn.interAreaDelay/2);
    upDelay = ceil(param.conn.interAreaDelay/2);
    
    m.add(synL, 'T', upDelay, 2);
    m.add(synR, 'T', upDelay, 2);
    m.add(fDnL, 'T', dnDelay, 2);
    m.add(fDnR, 'T', dnDelay, 2);
    
end

out = m.run('struct');


% get measures for hemispheric biases:

% stim-off state prob:
out.lower_bias = log(out.lower_L./out.lower_R);
out.upper_bias = log(out.upper_L./out.upper_R);




end


