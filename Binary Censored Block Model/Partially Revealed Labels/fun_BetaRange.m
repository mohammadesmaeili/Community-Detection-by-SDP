function [Beta_range] = fun_BetaRange(n, a, xi, M)
Beta = 0;
Erasure = n.^-Beta;
[p] = fun_parameters(n, a);
Eta0 = a.*(sqrt(1-xi)-sqrt(xi)).^2+Beta;

syms Beta
Erasure = n.^-Beta;
p = fun_parameters(n, a);
Eta = a.*(sqrt(1-xi)-sqrt(xi)).^2+Beta;

Beta_range = [];
for i=1:length(M)
    x = vpasolve(Eta==M(i),Beta);
    x = double(x);
    Beta_range = [Beta_range,x];
end