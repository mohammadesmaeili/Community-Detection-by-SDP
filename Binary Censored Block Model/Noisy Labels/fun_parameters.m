function [p, T1, T2]=fun_parameters(n, a, xi, alpha)
p = a*log(n)/n;
T1 = log((1-xi)/xi);
T2 = log((1-alpha)/alpha);