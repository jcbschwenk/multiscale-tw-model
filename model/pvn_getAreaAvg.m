function out = pvn_getAreaAvg(mf, g, param, flag)
% returns averages of all local activity in each area.
% Cx and Pulvinar are treated separately.
% input is:
% mf - output of meanfield model
% g - graph object
% param - parameter struct
% flag - (optional) can be 'Cx' or 'Pul' to return only Cx or Pul

if nargin < 4
   flag = ''; 
end

for iRun = 1:numel(mf)
    out(iRun).avg = [];
    out(iRun).label = {};
    for iLevel = 0:param.N+1
        idxCx = contains(g.Labels, num2str(iLevel)) & ~contains(g.Labels, 'Pul');
        idxPul = find(contains(g.Labels, num2str(iLevel)) & contains(g.Labels, 'Pul'));
        if ~isempty(idxPul) && ~strcmpi(flag, 'Pul')
            out(iRun).avg(end+1,:) = mean(mf(iRun).r(idxCx,:),1);
            out(iRun).label{end+1} = ['Cx' num2str(iLevel)];
        end
        if ~isempty(idxPul) && ~strcmpi(flag, 'Cx')
            out(iRun).avg(end+1,:) = mean(mf(iRun).r(idxPul,:),1);
            out(iRun).label{end+1} = ['Pul' num2str(iLevel)];
        end
    end
end

[out.t] = deal(mf.t);


end