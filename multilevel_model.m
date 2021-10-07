function [A, h] = multilevel_model(K,m,gamma,a)
% multilevel_model  Generate synthetic network from the trophic RDRG model
% 
% INPUTS
% - K Number of clusters
% - m Number of nodes per cluster
% - gamma  Decay parameter
% - a  Parameter for the additive noise
%
% OUTPUTS
% - A Adjacency matrix of the synthetic network
% - theta  True phase angles of nodes

n = m*K; % total number of nodes
h = sort(repmat(linspace(0,K-1,K),1,m)+(2*a*rand(1,n)-a)); % trophic levels
A = zeros(n,n); % preallocate adjacency matrix
f = zeros(n,n); % preallocate edge probability
F = zeros(n,n); % preallocate trophic incoherence

for i = 1: (n-1)
    for j = i+1 : n
        F(i,j) = (h(j)-h(i)-1)^2; % trophic incoherence of edge i->j
        F(j,i) = (h(i)-h(j)-1)^2; % trophic incoherence of edge j->i
        f(i,j) = 1/(1+exp(gamma* F(i,j))); % probability of edge i->j
        f(j,i) = 1/(1+exp(gamma* F(j,i))); % probability of edge j->i
        
        %generate edge i->j with probaiblity f(i,j)
        if rand <= f(i,j)
            A(i,j) = 1;
        else
            A(i,j) = 0;
        end
        
        %generate edge j->i with probaiblity f(j,i)
        if rand <= f(j,i)
            A(j,i) = 1;
        else
            A(j,i) = 0;
        end
    end
end        
end
