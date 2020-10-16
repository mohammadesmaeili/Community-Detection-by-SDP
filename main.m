clc
clear all
close all
format long
%% variables
n = 100;
rho = 0.5;
a = 4;
b = 1;
Beta = 0.5;
alpha = 1/(1+n^Beta);
MaxCounter = 1e8;
SaveCounter = 1e6;
%%
[p, q, T1, T2] = fun_parameters(n, a, b, alpha);
Gamma = sqrt(Beta^2+a*b*T1^2);
Eta = (a+b)/2+Beta/2-Gamma/T1+Beta/2/T1*log((Gamma+Beta)/(Gamma-Beta))
%%
Error = comm.ErrorRate;
counter = 0;
Check = 1;
error_dec = [];
norm_dec = [];
iter = 0;
mean_norm = 0;
plot_norm = [];
plot_error = [];
while Check
    [X, Y] = fun_generating_node_features(n,rho,alpha);
    G = fun_graph_generator(X,p,q);
    [x_hat, norm_diff] = SDP_Noisy_BSBM(G, X, Y, T1, T2);

    Er=Error(X,x_hat);
    Er(1)
    mean_norm = (iter*mean_norm+norm_diff)/(iter+1)
    counter = counter + n;
    iter = iter + 1;
    plot_norm = [plot_norm, mean_norm];
    plot_error = [plot_error, Er(1)];
    if length(plot_norm)>100
        plot_norm = plot_norm(2:end);
        plot_error = plot_error(2:end);
    end
    figure(1)
    plot(abs(plot_norm-mean_norm)/mean_norm)
    figure(2)
    plot(abs(plot_error-Er(1))/Er(1))
        
    if counter>SaveCounter
        Str=['Sim','_a10_',num2str(a*10),...
            '_b10_',num2str(b*10),...
            '_B_',num2str(Beta*10),...
            '_n_',num2str(n),...
            '_MaxBit_',num2str(Er(3))];
        save(Str)
        norm_dec = [norm_dec, mean_norm ]
        error_dec = [error_dec, Er(1)]
        if Er(3)>=MaxCounter
            Check = 0;
        end
        counter = 0;
    end
end
%%
Str=['Sim','_a10_',num2str(a*10),...
    '_b10_',num2str(b*10),...
    '_B_',num2str(Beta*10),...
    '_n_',num2str(n),...
    '_MaxBit_',num2str(Er(3))];
save(Str)
reset(Error)