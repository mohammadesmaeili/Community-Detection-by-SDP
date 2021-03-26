function [Beta_range] = fun_BetaRange(n, a, b, M)
Beta = 0;
alpha = 1./(1+n.^Beta);
[p, q, T1, T2] = fun_parameters(n, a, b, alpha);
Gamma = sqrt(Beta.^2+a.*b.*T1.^2);
Eta0 = (a+b)./2+Beta./2-Gamma./T1+Beta./2./T1.*log((Gamma+Beta)./(Gamma-Beta));

syms Beta
alpha = 1./(1+n.^Beta);
[p, q, T1, T2] = fun_parameters(n, a, b, alpha);
Gamma = sqrt(Beta.^2+a.*b.*T1.^2);
Eta = (a+b)./2+Beta./2-Gamma./T1+Beta./2./T1.*log((Gamma+Beta)./(Gamma-Beta));

Beta_range = [];
for i=1:length(M)
    x = vpasolve(Eta==M(i),Beta);
    x = double(x);
    Beta_range = [Beta_range,x];
end