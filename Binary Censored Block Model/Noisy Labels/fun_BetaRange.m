function [Beta_range] = fun_BetaRange(n, a, xi, M)
Beta = 0;
alpha = 1/(1+n^Beta);
[~, T1, ~] = fun_parameters(n, a, xi, alpha);
Gamma = sqrt(Beta.^2+4.*xi.*(1-xi).*a.^2.*T1.^2);
Eta0 = a-Gamma./T1+Beta./2./T1.*log((1-xi).*(Gamma+Beta)./xi./(Gamma-Beta));

syms Beta
alpha = 1/(1+n^Beta);
[~, T1, ~] = fun_parameters(n, a, xi, alpha);
Gamma = sqrt(Beta.^2+4.*xi.*(1-xi).*a.^2.*T1.^2);
Eta = a-Gamma./T1+Beta./2./T1.*log((1-xi).*(Gamma+Beta)./xi./(Gamma-Beta));

Beta_range = [];
for i=1:length(M)
    x = vpasolve(Eta==M(i),Beta);
    x = double(x);
    Beta_range = [Beta_range,x];
end