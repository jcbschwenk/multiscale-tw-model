function [ln, ppatch] = pvn_shplot(x,y,e,varargin)
keephold = ishold;
if iscolumn(e) || size(e,1) > size(e,2)
    e = e';
end
if iscolumn(y)
    y = y';
end
hold on
if size(e,1) == 2
    e_u = e(1,:);
    e_l = e(2,:);
else
    e_u = y + 0.5.*e; e_l = y - 0.5.*e;
end

ln = plot(x,y,varargin{:});
mainCol = 0.3*ln.Color + 0.7;

set(gcf,'renderer','painters');
ppatch = patch([x flip(x)],[e_l flip(e_u)],1,'FaceColor',mainCol,'edgecolor','none','facealpha',0.5);

uistack(ln,'top');
if ~keephold
    hold off
end

end