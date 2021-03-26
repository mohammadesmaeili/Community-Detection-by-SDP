function [x_hat, norm_diff]=SDP_Erasure_BCBM(A, x, Y)
n = length(x);
J=ones(n,n);
W = Y*Y';
cvx_begin sdp
    variable Z(n,n) symmetric
    maximize ( trace(A*Z) )
    Z >= 0
    diag(Z) == 1;
    trace( W*Z ) == (Y'*Y)^2;
cvx_end

norm_diff = norm(x*x'-Z)/n;
%% first decoder
x_hat = Z*x;
x_hat_1 = x_hat/n;
%% second decoder
x_hat = x_hat_1;
x_hat(x_hat_1>=0) = 1;
x_hat(x_hat_1<0) = -1;
end