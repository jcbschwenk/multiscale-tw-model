function [fig1, fig2] = pvn_singleCxStreamWaveDynamics()
% Runs simulation on single-stream, cortex-only version of the model.
% fig1: panel figure with mean-field responses, cortex average, eeg
% and wave fit; fig2: topography plots of mean phase gradients for BW and FW
% state.


[g, param, stim, prior] = pvn_mfmodel();
nTr = 50;
[stim, ~, isOnStim] = pvn_defaultStim('dcPulse_wnPrior', nTr, stim, prior);

%%
mf = g.run();
cx = pvn_getAreaAvg(mf, g, param, 'Cx');

eeg = pvn_eegProject(mf, g, param);
wav = pvn_fitPlaneEEG(eeg, 'Frequency', [7 13]);

%%
fig1 = figure;
pltTrl = 1;
cols = pvn_figCols('cx');

tiledlayout(6,1);

lblargs = {'FontWeight', 'Bold', 'FontName', 'Arial', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top'};

% Stimulus trace: 
nexttile
hold on
plot(stim.t, stim.I(1,:,pltTrl), 'k', 'LineWidth', 1.5);
txtpos(1) = min(xlim) + diff(xlim)*0.05;
text(txtpos(1), max(ylim)*0.9, 'STIM', lblargs{:})
axis off

% SGX response traces: 
nexttile
hold on
for i = 1:param.N
    plot(mf(pltTrl).t, mf(pltTrl).r(g.getNodeIdx(['SGX' num2str(i)]),:), 'Color', cols{i}, 'LineWidth', 1.5);
end
text(txtpos(1), max(ylim)*0.9, 'SGX', lblargs{:})
axis off
% patch([isOnStim{1} flip(isOnStim{1})], [0 1; 0 1; 1 0; 1 0] * [min(ylim) max(ylim)]', [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% IGX (IGIB) response traces: 
nexttile
hold on
for i = 1:param.N
    plot(mf(pltTrl).t, mf(pltTrl).r(g.getNodeIdx(['IGX' num2str(i)]),:), 'Color', cols{i}, 'LineWidth', 1.5);
end
text(txtpos(1), max(ylim)*0.9, 'IGX', lblargs{:})
axis off

% Heatmap of cortex averages:
nexttile
imagesc(cx(pltTrl).t, 1:param.N, cx(pltTrl).avg)
xline(0, 'w:', 'LineWidth', 2)
xlabel('Time [ms]')
yticks([1 2 3])
yticklabels(cx(pltTrl).label)
text(txtpos(1), max(ylim)*-0.1, 'Cx Avg', lblargs{:})
colormap(brewermap(64, '*RdYlBu'))
cb = colorbar;
cb.Ticks = cb.Ticks([1 end]);


% Heatmap of single-trial EEG (midline):
nexttile
imagesc(eeg.t, 1:numel(eeg.midlineIdx), eeg.eeg(eeg.midlineIdx,:,pltTrl))
xline(0, 'w:', 'LineWidth', 2)
yticks(1:numel(eeg.midlineIdx))
yticklabels(eeg.label(eeg.midlineIdx))

% Trace of single-trial wave-state, with marked epochs:
nexttile
plot(wav.t, wav.wavDir(:,pltTrl), 'k', 'LineWidth', 1.5);
hold on
wdirs = {'fw' 'bw'};
stateEpochs = cellfun(@(x) bwlabel(wav.sig(:,pltTrl) & wav.(x)(:,pltTrl)), wdirs, 'UniformOutput', false);
stateEpochs = cellfun(@(x) arrayfun(@(y) [find(x == y,1,'first') find(x == y,1,'last')], unique(x(x > 0)), 'UniformOutput', false), stateEpochs, 'UniformOutput', false);
for iwd = 1:2
    for iEpoch = 1:numel(stateEpochs{iwd})
        xx = wav.t(stateEpochs{iwd}{iEpoch});
        if diff(xx) < 0.2
           continue 
        end
        patch([xx flip(xx)], [-1 -1 1 1].*pi, pvn_figCols(wdirs{iwd}), 'EdgeColor', 'none');
        alpha(0.2)
    end
end
xlabel('Time [ms]')
yline(pi/2, 'k:')
yline(-pi/2, 'k:')
ylim([-pi pi])
yticks([-pi/2 pi/2])
yticklabels({'FW' 'BW'})


set(gcf, 'Color', 'w');

%%
fig2 = figure;

s1 = subplot(1,2,1);
pvn_plotEEGTopo(eeg, 'TimeWindow', [-2 -1], 'Parameter', 'Phase', 'PlotType', '2D', 'PhaseReference', 'Iz', 'FigHandle', s1);
title('BW')
s2 = subplot(1,2,2);
pvn_plotEEGTopo(eeg, 'TimeWindow', [0 1], 'Parameter', 'Phase', 'PlotType', '2D', 'PhaseReference', 'Iz', 'FigHandle', s2)
title('FW')
cb = colorbar('Location', 'EastOutside');
cb.Position = cb.Position + [0.1 0 0 0]; % shift colorbar to get original subplot scale back


end