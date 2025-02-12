function r = pvn_ptBiSerialCorr(dich, cont)
% computes the point-biserial correlation between dichotomous variabl dich
% and continuous variable cont

assert(numel(unique(dich)) == 2, 'First input variable must be dichotomous')

assert(numel(dich) == numel(cont), 'Dimensions of input variables mismatch')

vals = unique(dich);

for i = 1:2
   m(i) = mean(cont(dich == vals(i)));
   n(i) = sum(dich == vals(i));
end
s = std(cont);
nTot = numel(cont);

r = (diff(m)/s) * sqrt(prod(n)/nTot^2);

end