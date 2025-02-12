function out = pvn_fitPlaneEEG(eeg, varargin)
% Plane fit for TW classification. 
% Based on tw_fitPlaneEEG, modified for pvn.
%
% input 'eeg' is output of pvn_eegProject;
%
% Optional Input (as named arg list):
%
% 'Frequency'           temporal freq limits [fmin fmax] used for bp filter
% 'ROI'                 defines limits of electrode region to be used for the fit 
%                       (see fieldtrip layout EEG1005)
% 'WindowSize'          temporal window size (ms) used to average relative phase 
%                       (leave empty if single-point)
% 'MaxCycles'           max. spatial frequency, defined as number of cycles
%                       within the covered electrode space, default 1 cycle
% 'NumStepsSpatFreq'    number of steps for spatial freq [0 MaxCycles],
%                       default 30
% 'NumStepsWaveDir'     number of steps for wave direction [-pi pi],
%                       default 60
% 'RandShuffleIter'     number of iterations for random
%                       shuffling to estimate the chance level of
%                       rcc_square.
% 'RandShuffleNTrials'  number of trials to perform random shuffling on. If
%                       empty, all trials 
% 'DirClassWindowSize'  window size for classification of wave direction (in radians)
% 'AdjustDirAxis'       if true, wave dir is classified relative to median
%                       of fitted wave directions (after mirroring along M-L); 
%                       if false, relative to [+/- pi/2]. Use if tilted
%                       axis is expected (e.g. hemispheric fit) and only if
%                       median estimation is expected to be accurate.
% 'RandShuffleTimeStep' time step option for shuffling, in ms. 

p = inputParser();
p.addParameter('Frequency', [7 13]); 
p.addParameter('ROI', {[-0.2 0.2] [-Inf 0.25]});
p.addParameter('WindowSize', []); 
p.addParameter('MaxCycles', 1); 
p.addParameter('NumStepsSpatFreq', 30); 
p.addParameter('NumStepsWaveDir', 60); 
p.addParameter('RandShuffleIter', 10); 
p.addParameter('RandShuffleNTrials', []);
p.addParameter('DirClassWindowSize', 0.5); 
p.addParameter('AdjustDirAxis', false); 
p.addParameter('RandShuffleTimeStep', 200); % in ms

p.parse(varargin{:});

%%
raw = eeg.eeg;
nTr = size(raw,3);
nTm = size(raw,2);
sr = round(1/diff(eeg.t(1:2)));

% convert moving-avg windows size to samples:
if ~isempty(p.Results.WindowSize)
    movMeanWinSize = (p.Results.WindowSize/1e3) .* sr;
else
    movMeanWinSize = [];
end

if ~isempty(p.Results.RandShuffleTimeStep)
    shStep = round((p.Results.RandShuffleTimeStep/1e3) .* sr);
    
    shTm = shStep:shStep:nTm;
else
    shTm = 1:nTm;
end

if isempty(p.Results.RandShuffleNTrials)
    shuffleTrials = 1:nTr;
else
    shuffleTrials = sort(randsample(1:size(raw,3),...
        min(nTr, p.Results.RandShuffleNTrials)));
end

%%

% Get the electrode positions that we want to include in the fit:
[pos, lbl, idxInData] = getElecPos(eeg.label, 'EEG1005', p.Results.ROI);


% Get the indices for the midline electrodes (not used in the fit):
[midLine, midLineLabel, midLineIdxInData] = getMidLine(lbl, idxInData);

% Filter the data and extract phases:
[phi, pow, eegFilt] = getBandpassPhase(raw(idxInData,:,:), sr, p.Results.Frequency);

% Get the phase planes predicted by each fit:
[phiPred, wvdir, a, b, xi] = getPlaneFits(pos,...
    p.Results.NumStepsSpatFreq, p.Results.NumStepsSpatFreq, p.Results.MaxCycles);
% -> first dimension for each is the number of different fits

% Evaluate the fits:
[id, rcc_sq, rcc_sq_rand, phi_out] = evalPlaneFits(phi, phiPred, movMeanWinSize, p.Results.RandShuffleIter, shuffleTrials, shTm);
 
[fw, bw] = classifyDirection(wvdir(id), p.Results.DirClassWindowSize, p.Results.AdjustDirAxis);

%%
out.t = eeg.t;
out.fw = fw;
out.bw = bw;

out.lbl = lbl;
out.pos = pos;
out.xi = xi(id);
out.a = a(id);
out.b = b(id);
out.wavDir = wvdir(id);

out.rcc_sq = rcc_sq;
out.rcc_sq_rand = rcc_sq_rand(:,shuffleTrials,:);
out.rcc_sq_rand = out.rcc_sq_rand(shTm, :,:);
out.rand_trials = shuffleTrials;
out.rand_t = out.t(shTm);
out.rcc_thresh = prctile(out.rcc_sq_rand(:), 95);
out.sig = out.rcc_sq > out.rcc_thresh;

out.pFitFW = mean(out.sig & out.fw,2);
out.pFitBW = mean(out.sig & out.bw,2);
out.pFitFW_CI = bootStrapProbCI(out.fw, out.sig, 100);
out.pFitBW_CI = bootStrapProbCI(out.bw, out.sig, 100);

out.midLineIdx = midLine;
out.midLineIdxInData = midLineIdxInData;
out.midLineLabel = midLineLabel;
out.pow = pow;
out.phi = phi_out;
out.label  = lbl;

%%
end
%%
function [pos, lbl, idxInData] = getElecPos(labels, layout, boxLims)


lay = ft_prepare_layout(struct('layout', layout));
box = @(xy) xy(:,1) >= boxLims{1}(1) & xy(:,1) <= boxLims{1}(2) & xy(:,2) >= boxLims{2}(1) & xy(:,2) <= boxLims{2}(2);

inclFromLay = find(box(lay.pos));
idxInData = find(cellfun(@(x) ismember(x, lay.label(inclFromLay)), labels));

idxInLay = nan(1,numel(idxInData));
for iElec = 1:numel(idxInData)
    idxInLay(iElec) = find(strcmpi(lay.label, labels{idxInData(iElec)}));
end

pos = lay.pos(idxInLay,:);
pos = pos - mean(pos); % center on zero/zero
lbl = labels(idxInData);
end
%%
function [id, rcc_sq, rcc_sq_rand, phi_out] = evalPlaneFits(phi, phiPred, movMeanWinSize, nIter, shuffleTrials, shTm)
% compares predicted phases in phiPred to actual phases in phi, after
% applying moving-mean of window size movMeanWinSize (in samples, not ms);
% phases in phi must be relative (e.g. to mean phase across electrode) if
% moving mean is requested
% Returned id is the index of the best fit for each time-point, with the same
% [nTime x nTrial] dimensions as input

%%
nTr = size(phi,3);
nTm = size(phi,2);

%%
if nargin < 3 || isempty(movMeanWinSize)
    movAvgIdx = @(t) t;
else
    movAvgIdx = @(t) (1:nTm) > t-movMeanWinSize/2 & (1:nTm) < t+movMeanWinSize/2;
end
if nargin < 4 || isempty(nIter)
    nIter = 0;
end

%% Random shuffling:
if size(phi,1) >= 10
    for iter = 1:nIter
        phi_rand(:,:,:,iter) = phi(randperm(size(phi,1)),:,:);
    end
else
    assert(nIter < factorial(size(phi,1)),...
        'Number of iterations is greater than number of possible permutations')
    rperm = perms(1:size(phi,1));
    rperm = rperm(randsample(1:size(rperm,1), nIter),:);
    
    for iter = 1:nIter
        phi_rand(:,:,:,iter) = phi(rperm(iter,:),:,:);
    end
    
end

%%
id = nan(nTm,nTr);
rcc_sq = nan(nTm,nTr);
rcc_sq_rand = nan(nTm,nTr,nIter);
phi_out = nan(size(phi,1), nTm, nTr);

wb = waitbar(0, {sprintf('Wave Fit: Trial %d/%d', 1,nTr), '~~~ >> ??? << ~~~'});

for itrial = 1:nTr
    wbstr = {sprintf('Wave Fit: Trial %d/%d', itrial,nTr), '~~~ >> ??? << ~~~'};
%     if ~exist('wb') || ~isvalid(wb)
%         wb = waitbar(itrial/nTr, wbstr);
%     else
        waitbar(itrial/nTr,wb, wbstr);
%     end
    for iT = 1:nTm
        aphi = squeeze(circ_mean(phi(:,movAvgIdx(iT),itrial), [], 2));
        [id(iT,itrial), rcc_sq(iT,itrial)] = singleFitEval(aphi, phiPred);
        phi_out(:,iT,itrial) = aphi;
    end
    
    if ismember(itrial, shuffleTrials)
        for iT = shTm
            for iter = 1:nIter
                rphi = squeeze(circ_mean(phi_rand(:,movAvgIdx(iT),itrial,iter), [], 2));
                [~, rcc_sq_rand(iT,itrial,iter)] = singleFitEval(rphi, phiPred);
            end
        end
    end
    
end

delete(wb);

phi_out = phi_out - circ_mean(phi_out, [], 1); % recenter on zero after moving mean

end
%%
function [id, rcc_sq] = singleFitEval(aphi, phiPred)
gof = sqrt(mean(cos(phiPred'-aphi)).^2 + mean(sin(phiPred'-aphi)).^2); % mean vector length of residuals (-> ie offset is ignored)
[~,id] = max(gof); % get best fit ID
bestfit = phiPred(id,:)'; % get predicted phases for best fit
rcc_sq = sum(sin(aphi-circ_mean(aphi)) .* sin(bestfit-circ_mean(bestfit)))...
    ./sqrt(...
    sum(sin(aphi-circ_mean(aphi)).^2) .* sum(sin(bestfit-circ_mean(bestfit)).^2)); % circular correlation for best fit
rcc_sq = rcc_sq^2;
end
%%
function [phi, pow, eegFilt] = getBandpassPhase(raw, sr, freq)
% bandpass filter raw data and return phases
% raw is raw eeg signal as [nChan x nTime x nTrials]
% freq is band as [min max] in Hz

%%
nTr = size(raw,3);
nTm = size(raw,2);
nElec = size(raw,1);

%% Filter the data
bpFilt = designfilt('bandpassfir', 'StopbandFrequency1', 0.85*freq(1), 'PassbandFrequency1', freq(1),...
    'PassbandFrequency2', freq(2), 'StopbandFrequency2', 1.25*freq(2),...
    'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', sr);


eegFilt = nan(nElec, nTm, nTr);

try
    for itrial = 1:nTr
        for iElec = 1:nElec
            eegFilt(iElec,:,itrial) = filtfilt(bpFilt, raw(iElec,:,itrial));
        end
    end
catch errmsg
    if strcmpi(errmsg.identifier, 'signal:filtfilt:InvalidDimensionsDataShortForFiltOrder')
        rawPad = cat(2,raw,raw,raw);
        eegFilt = nan(nElec, 3*nTm, nTr);
        for itrial = 1:nTr
            for iElec = 1:nElec
                eegFilt(iElec,:,itrial) = filtfilt(bpFilt, rawPad(iElec,:,itrial));
            end
        end
        eegFilt = eegFilt(:,nTm+1:2*nTm,:);
    else
        throw(errmsg)
    end
end


%% extract phase:
for itrial = 1:nTr
    h = hilbert(eegFilt(:,:,itrial)');
    phi(:,:,itrial) = angle(h)';
    pow(:,itrial) = mean(abs(h),2);
end
phi = phi - circ_mean(phi, [], 1); % subtract mean across elecs per timepoint
% -> now we can average over time (within short windows)


end
%%
function [midLine, midLineLabel, midLineIdxInData] = getMidLine(lbl, idxInData)
% get indices for midline electrodes for a given set of electrode labels
% passed as cell array

midLineLabel = {'Iz' 'Oz', 'POz', 'Pz', 'CPz', 'Cz', 'FCz', 'Fz'};

for i = 1:numel(midLineLabel)
    idx = find(strcmpi(lbl, midLineLabel{i}));
    if ~isempty(idx)
        midLine(i) = idx;
        midLineIdxInData(i) = idxInData(idx);
    else
        midLine(i) = nan;
        midLineIdxInData(i) = nan;
    end
end
midLine = midLine(~isnan(midLine));
midLineLabel = midLineLabel(~isnan(midLine));

end
%%
function [phiPred, wvdir, a, b, xi] = getPlaneFits(pos, numStepsWaveDir, numStepsSpatFreq, numCyclesMax)
% get predicted phases for all plane fits. NumStepsWaveDir / SpatFreq are
% the resolution in radial and circular dimensions, numCyclesMax is the
% maximum radius (i.e. max spatial freq.), as number of cycles across the
% covered electrode space

deltapos = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2);
d = max(deltapos(:)); % this is the max distance covered by the electrode space

xiMax = (numCyclesMax*2*pi)/d;

xiStep = xiMax/numStepsSpatFreq;
thetaStep = (2*pi)/numStepsWaveDir;

theta = -pi:thetaStep:pi-thetaStep; % circular parameter (-> wave direction)
a = (xiStep:xiStep:xiMax)' .* sin(theta);
b = (xiStep:xiStep:xiMax)' .* cos(theta);

a = a(:);
b = b(:);

xi = sqrt(a.^2 + b.^2); % radial parameter / spatial freq
wvdir = atan2(b,a);
phiPred = a.*pos(:,1)' + b.*pos(:,2)'; % predicted phases
end
%%
function [fw, bw] = classifyDirection(wavDir, winSize, useMedian)

if useMedian
    % mirror along the M-L axis to get a shared histogram: (this assumes
    % that there is a single axis of propagation, that is closer to the A-P
    % axis! That's true for the model, but possibly not for real data...)
    
    wdDist = [wavDir(wavDir > 0); pi-abs(wavDir(wavDir < 0))];
    
    md = median(wdDist);
    
    fwMode = -pi+md;
    bwMode = md;
    
else
    fwMode = -pi/2;
    bwMode = pi/2;
end

fw = abs(wavDir-fwMode) < winSize;
bw = abs(wavDir-bwMode) < winSize;

end
%%
function ci = bootStrapProbCI(p, sig, nIter)
% estimates the two-sided 95% quantile for the state probabilites by
% bootstrapping (nIter samples with replacement).
% p is the trial-dim output for the raw classified state (fw/bw), sig is
% trial-dim significance based on rcc_sq threshold

N = size(p,2);

out = nan(size(p,1), nIter);

for iter = 1:nIter
   smpl = randsample(1:N, N, true); 
   out(:,iter) = mean(p(:,smpl) & sig(:,smpl),2);
end

ci = quantile(out, [0.025 0.975], 2);

end