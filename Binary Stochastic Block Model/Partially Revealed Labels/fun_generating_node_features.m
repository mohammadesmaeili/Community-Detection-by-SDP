function [X, Y] = fun_generating_node_features(n,rho, Erasure)
K=ceil(n*rho);
X=[ones(K,1);-ones(n-K,1)];
permx = randperm(n);
X=X(permx);
R = binornd(1,1-Erasure,n,1);
Y = R.*X;
end

