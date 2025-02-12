function out = pvn_getWTSpec(in, g, varargin)

p = inputParser();
p.addParameter('TimeWindow', [-Inf Inf]);
p.addParameter('Nodes', g.Labels); % can be labels or node types or idx
p.addParameter('FreqBandPeak', [0 Inf]);
p.addParameter('Field', 'r');

p.parse(varargin{:});

fld = p.Results.Field;
f = 1:0.1:100;
spec = nan(size(in(1).(fld),1), numel(f), numel(in));

for iTr = 1:numel(in)
    
    tIdx = in(iTr).t >= p.Results.TimeWindow(1) & in(iTr).t < p.Results.TimeWindow(2);
    
    for iNode = 1:size(in(iTr).(fld),1)
        [wt, wtFreq] = cwt(in(iTr).(fld)(iNode,:), 1e3, 'FrequencyLimits', [1 100]);
        wt = wt(end:-1:1,:);
        spec(iNode,:,iTr) = interp1(flip(wtFreq), mean(abs(wt(:,tIdx)),2), f);
    end
end

out.f = f;
bandlim = @(x) x >= p.Results.FreqBandPeak(1) & x < p.Results.FreqBandPeak(2);

for iType = 1:numel(p.Results.Nodes)
    if ~isnumeric(p.Results.Nodes)
        thisType = cellfun(@(x) contains(x, p.Results.Nodes{iType}), g.Labels);
    else
        thisType = iType;
    end
    
    for iTr = 1:numel(in)
        out.spec(iType,:,iTr) = mean(spec(thisType,:,iTr),1);
        [pks, locs] = findpeaks(out.spec(iType,:,iTr),out.f);
        if ~any(bandlim(locs))
            out.pk(iType,iTr) = nan;
            continue
        end
        pks = pks(bandlim(locs));
        locs = locs(bandlim(locs));
        [~, idx] = max(pks);
        out.pk(iType,iTr) = locs(idx);
    end
end

out.label = p.Results.Nodes;

end