function fig = pvn_dualStreamWaveDynamics()
% Runs simulation on dual-stream, cortex-only version of the model.

[g, param, stim, prior] = pvn_mfmodel('dual');

nTr = 50;

[~, ~, isOnStim] = pvn_defaultStim('dcPulseLeft_sharedWnPrior_sharedWnStim', nTr, stim, prior);

%%
mf = g.run();
eeg = pvn_eegProject(mf, g, param);
wavLeft = pvn_fitPlaneEEG(eeg, 'Frequency', [7 13], 'ROI',  {[-0.3 -0.05] [-Inf 0.25]});
wavRight = pvn_fitPlaneEEG(eeg, 'Frequency', [7 13], 'ROI',  {[0.05 0.3] [-Inf 0.25]});

%%

%%
fig = figure;

subplot(1,2,1)
ln_left = pvn_shplot(wavLeft.t, wavLeft.pFitFW, wavLeft.pFitFW_CI, 'Color', pvn_figCols('fw') + [0.4 0.4 0.3]);
hold on
ln_right = pvn_shplot(wavRight.t, wavRight.pFitFW, wavRight.pFitFW_CI, 'Color', pvn_figCols('fw') - [0.2 0.2 0.4]);
ylim([-0.1 1.1])
yline(0, 'k:');
yline(1, 'k:');
patch([isOnStim{1} flip(isOnStim{1})]./1e3, [-0.1 -0.1 1.1 1.1], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

xlabel('Time [secs]')
ylabel('P(FW)')
legend([ln_left ln_right], {'L' 'R'})
legend('boxoff')

s2 = subplot(1,2,2);
pvn_plotEEGTopo(eeg,...
    'TimeWindow', [0 0.1],...
    'Parameter', 'RMS',...
    'PlotType', '2D',...
    'ROIMask', {[-Inf Inf] [-Inf Inf]},...
    'FigHandle', s2);

set(gcf, 'Color', 'w')

%%
end