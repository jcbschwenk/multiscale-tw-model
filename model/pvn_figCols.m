function out = pvn_figCols(req)
% defines standard colors for different conditions between figures


out.fw = [0.2 0.2 0.7];
out.bw = [0.85 0.1 0.15];

out.cx = {[0.29 0.58 0.29] [0.65 0.58 0.83] [0.9 0.7 0.3]};

out.fwL = [0.45 0.45 1];
out.fwR = [0 0 0.4];

out.bwL = [0.9 0.1 0.15];
out.bwR = [0.65 0.15 0.2];



if nargin
   out = out.(req); 
end

end