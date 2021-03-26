function [Beta_range] = fun_BetaRange(n, a, b, M)
Beta = 0;
[p, q] = fun_parameters(n, a, b);
Eta0 = 0.5.*(sqrt(a)-sqrt(b)).^2 + Beta;

syms Beta
[p, q] = fun_parameters(n, a, b);
Eta = 0.5.*(sqrt(a)-sqrt(b)).^2 + Beta;

Beta_range = [];
for i=1:length(M)
    x = vpasolve(Eta==M(i),Beta);
    x = double(x);
    Beta_range = [Beta_range,x];
end