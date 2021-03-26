function [x_hat, norm_diff]=SDP_Noisy_BSBM(A, x, Y, T1, T2)
n = length(x);
J=ones(n,n);
cvx_begin sdp
    variable Z(n,n) symmetric
    variable X(n,1)
    maximize ( T1*trace(A*Z)+2*T2*X'*Y )
    [1, X' ; ...
        X, Z] >= 0
    diag(Z) == 1;
    trace( J*Z ) == 0;
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