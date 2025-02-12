function out = pvn_eegProject(in, g, param, src_areas)

isdual = any(contains(param.NodeLabels, '_R')); % dual hemisphere model, i.e. we have independent sources for each side

%%
if nargin < 4
   src_areas = param.eeg.source_labels; 
end
assert(numel(src_areas) == param.N || numel(src_areas) == param.N*2, 'Wrong number of input areas...');

srcPos = cellfun(@(x) param.eeg.source_pos(strcmpi(param.eeg.source_labels, x)), src_areas);
srcPos = cat(1, srcPos{:});

%%
assert(param.N < 8, 'Automatic area assignment assumes single digit labels.')
whichArea = cellfun(@(x) str2double(x(end)), g.Labels);  % get area index for each node

inclIdx = whichArea >= 1 & whichArea <= param.N; % exclude input-layers
inclIdx = find(inclIdx & ~contains(g.Labels, 'Pul')); % exclude Pulvinar
if isdual
    whichHemi = contains(g.Labels, '_R');
    srcPos = srcPos(whichArea(inclIdx)+whichHemi(inclIdx).*size(srcPos,1)/2, :); % repmat source positions
else
    srcPos = srcPos([whichArea(inclIdx) whichArea(inclIdx)+max(whichArea(inclIdx))], :);
end
%%
for iTrial = 1:numel(in)
    [r{iTrial}, t{iTrial}] = processSource(in(iTrial), inclIdx, param); % low-pass filter and resample
end
r = cellfun(@(x) x./rms(x(:)), r, 'UniformOutput', false); % normalize

if ~isdual
    % single-hemisphere model: 
    % we split each source between both hemispheres
    r = cellfun(@(x) 0.5.*[x; x], r, 'UniformOutput', false); 
end

for iTr = 1:numel(r)
    ns{iTr} = make_some_noise(param.eeg.num_noise_src, numel(t{iTr}));
    nsSNR = rand(param.eeg.num_noise_src, 1) .* diff(param.eeg.snr_range) + param.eeg.snr_range(1);
    ns{iTr} = (ns{iTr}./rms(ns{iTr},2)) ./ nsSNR;
end

nPos = size(param.eeg.ns_pos, 1);
nspos = randi(nPos, param.eeg.num_noise_src, numel(r));


src = cellfun(@(x,y) cat(1,x, y), r, ns, 'UniformOutput', false);

%%

for iTr = 1:numel(src)
    cfg = [];
    
    thisPos = cat(1, srcPos, param.eeg.ns_pos(nspos(:,iTr), :));
    
    cfg.sourcemodel.pos = thisPos;
    cfg.sourcemodel.mom = (thisPos/norm(thisPos))';
    cfg.sourcemodel.signal  = src(iTr);
    cfg.sourcemodel.time = t(iTrial);
    cfg.elec = param.eeg.elec;
    cfg.headmodel = param.eeg.headmodel;
    data = ft_dipolesimulation(cfg);
    
    
    out.eeg(:,:,iTr) = cat(3,data.trial{:});
end

out.label = data.label;
out.midlineIdx = cellfun(@(x) find(strcmpi(x, out.label)), param.eeg.midline);
out.t = t{1};



end
%%
function [out, t] = processSource(in, inclIdx, param)
    for iChan = 1:numel(inclIdx)
        tmp = filtfilt(param.eeg.bpFilt.b, param.eeg.bpFilt.a, in.r(inclIdx(iChan),:));
        [out(iChan,:), t] = resample(tmp, in.t./1e3, param.eeg.sr);
    end
end
%%
function out = make_some_noise(nSrc, nTm)
out = rand(nSrc,nTm);
for iSrc = 1:nSrc
    [wt, f] = cwt(out(iSrc,:));
    wt = wt./abs(wt);
    out(iSrc,:) = icwt(wt./sqrt(f));
end
end