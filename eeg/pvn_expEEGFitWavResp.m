% legacy function - not used for the final paper

function [ft, pltFt] = pvn_expEEGFitWavResp(varargin)
% fits a [rect x exp decay] model to the log-ratio of fw/bw wave
% probabilities (across trials). Assumes time windows based on the paradigm
% in Pang et al. 2020, used for real and simulated data.
% input is (t, logratio, isExp) or (wav, isExp) where wav is the output of
% pvn_fitPlaneEEG;
% if isExp is true, time window for fit anchor is shifted by 0.2 sec to
% account for latencies.

if isstruct(varargin{1}) && nargin == 2
    wav = varargin{1};
    t = wav.t;
    logp = log(wav.pFitFW./wav.pFitBW);
    isExp = varargin{2};
elseif ~isstruct(varargin{1}) && nargin == 3
    t = varargin{1};
    logp = varargin{2};
    isExp = varargin{3};
end

ftStartPtWin = [0.3 0.8]; % exp fit is anchored to the mean of this time window
ftStartPtWin = ftStartPtWin + isExp*0.2;
ftEndPtWin = [4 5]; % exp fit is anchored to the mean of this time window
ftEndPtWin = ftEndPtWin + isExp*0.2;
ftFitWin = [0.5 5]; % exp fit time window

trng = @(lm) t >= lm(1) & t < lm(2);

anchorStart = mean(logp(trng(ftStartPtWin) & ~isinf(logp')));
anchorEnd = mean(logp(trng(ftEndPtWin) & ~isinf(logp')));

bl = mean(logp((trng([-Inf 0]) | trng([6 Inf])) & ~isinf(logp')));

ftType = fittype('(hi-lo)*(1-rt)^x + lo', 'Independent', 'x', 'Dependent', 'y' );
ftOpts = fitoptions( 'Method', 'NonlinearLeastSquares');
ftOpts.Lower = [anchorStart anchorEnd -0.5];
ftOpts.Upper = [anchorStart anchorEnd 0.5];

t = t(~isinf(logp));
trng = @(lm) t >= lm(1) & t < lm(2);
logp = logp(~isinf(logp));

x = t(trng([0.5 5]))';
y = logp(trng(ftFitWin));

if numel(x) > 4
    ftRez = fit(x, y, ftType, ftOpts);
    pltFt = trng(ftFitWin).*ftRez(t)' + ~trng(ftFitWin).*bl; % create a version of the fit to plot
    
else
    [ftRez.hi, ftRez.lo, ftRez.rt] = deal(nan);
    pltFt = nan(size(t));
    
end

    
    ft.hi = ftRez.hi;
    ft.lo = ftRez.lo;
    ft.rt = ftRez.rt;
    ft.bl = bl;
    

end