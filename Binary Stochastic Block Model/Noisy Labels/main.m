clc
clear all
close all
format long
%% variables
n1 = 100:100:300;
rho = 0.5;
a = 3;
b = 1;
Beta_range = 0;
%%
Error = zeros(length(n1),1);
Er = zeros(length(n1),1);
MeanNorm = zeros(length(n1),1);
Check = 1;
iter = 0;
Eta_plot = zeros(length(Beta_range),1);
for i=1:length(n1)
    n = n1(i);
    Beta = Beta_range;
    alpha = 1/(1+n^Beta);
    [p, q, T1, T2] = fun_parameters(n, a, b, alpha);
    Gamma = sqrt(Beta^2+a*b*T1^2);
    Eta = (a+b)/2+Beta/2-Gamma/T1+Beta/2/T1*log((Gamma+Beta)/(Gamma-Beta));
    Eta_plot(i) = Eta;
end
Eta_plot
while Check
    iter = iter + 1; 
    parfor i=1:length(n1)
        %%
        n = n1(i);
        Beta = Beta_range;
        alpha = 1/(1+n^Beta);
        [p, q, T1, T2] = fun_parameters(n, a, b, alpha);
        Gamma = sqrt(Beta^2+a*b*T1^2);
        Eta = (a+b)/2+Beta/2-Gamma/T1+Beta/2/T1*log((Gamma+Beta)/(Gamma-Beta));
        %%
        [X, Y] = fun_generating_node_features(n,rho,alpha);
        G = fun_graph_generator(X,p,q);
        [x_hat, norm_diff] = SDP_Noisy_BSBM(G, X, Y, T1, T2);
        if X'*x_hat~=n
            Error(i) = Error(i)+1;
            Er(i) = Er(i) + (n-X'*x_hat)
        end
        MeanNorm(i) = ((iter-1)*MeanNorm(i)+norm_diff)/iter
        iter
    end
    figure(1)
    plot(Eta_plot, Error/iter, '-o')
    grid on
    figure(2)
    plot(Eta_plot, MeanNorm, '-*')
    grid on
    figure(3)
    semilogy(Eta_plot, Er'./iter./n1, 'o')
    grid on
    pause(0.1)
    if iter>=1000 && mod(iter,1000)==0
        Str=['Sim_new','_a_',num2str(a*10),...
            '_b_',num2str(b*10),...
            '_Beta_',num2str(Beta_range*10),...
            '_iter_',num2str(iter)];
        save(Str)
    end
end