function [p, q, T1, T2]=fun_parameters(n, a, b, alpha)
p = a*log(n)/n;
q = b*log(n)/n;
T1 = log(p*(1-q)/q/(1-p));
T2 = log((1-alpha)/alpha);