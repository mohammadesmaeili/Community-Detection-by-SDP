function [X, Y] = fun_generating_node_features(n,rho, alpha)
K=ceil(n*rho);
X=[ones(K,1);-ones(n-K,1)];
permx = randperm(n);
X=X(permx);
R = 2*binornd(1,1-alpha,n,1)-1;
Y = R.*X;
end

