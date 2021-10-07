function [W, h] = multilevel_model(k,n,gamma, a)
N = n*k;
h = sort(repmat(linspace(0,k-1,k),1,n)+(2*a*rand(1,N)-a));
%create empty adjacency matrix
W = zeros(N,N);

%create empty matrix for probability matrix
f = zeros(N,N); %edge probability
F = zeros(N,N); %incoherence

for i = 1: (N-1)
    for j = i+1 : N
        F(i,j) = (h(j)-h(i)-1)^2; %incoherence of edge i->j
        F(j,i) = (h(i)-h(j)-1)^2; %incoherence of edge j->i
        f(i,j) = 1/(1+exp(gamma* F(i,j))); %probability of edge i->j
        f(j,i) = 1/(1+exp(gamma* F(j,i))); %probability of edge j->i
        if rand <= f(i,j)
            W(i,j) = 1;
        else
            W(i,j) = 0;
        end
        
        if rand <= f(j,i)
            W(j,i) = 1;
        else
            W(j,i) = 0;
        end
    end
end        
end