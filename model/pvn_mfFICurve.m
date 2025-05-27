function y = pvn_mfFICurve(I, lm, b, c)
% FI curve for meanfield dynamics
y = (lm.*I - b)./(1- exp(-c.*(lm.*I-b)));

limIdx = (lm.*I - b) < 1e-8;
y(limIdx) = 1./c(limIdx);
end