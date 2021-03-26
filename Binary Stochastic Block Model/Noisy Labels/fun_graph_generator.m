function [G] = fun_graph_generator(X,p,q)
n = length(X);
G = zeros(n,n);
for i=1:length(X)
    for j=i+1:length(X)
        if X(i)==X(j)
            G(i,j)=binornd(1,p);
            G(j,i)=G(i,j);
        else
            G(i,j)=binornd(1,q);
            G(j,i)=G(i,j);
        end
    end
end
end

