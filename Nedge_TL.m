function [n_edge] = Nedge_TL(h, gamma)
% Nedge_TL   Calculate expected number of edges of the Trophic RDRG
%   
% - h       Estimated trophic levels of the nodes
% - gamma   Decay parameter for the Trophic RDRG model
%
% OUTPUTS
% - n_edge Expected number of edges

n_nodes = length(h);
n_edge = 0; %initialize number of edges
Isq = zeros(n_nodes, n_nodes); % preallocated trophic incoherence
z = zeros(n_nodes, n_nodes); % preallocate normalization factors
p = zeros(n_nodes, n_nodes); % preallocate probabilities of edges
for i = 1 : n_nodes - 1
    for j = i + 1 : n_nodes %loop through all pairs of nodes
        %calculate probability of edge
        Isq(i, j) = (h(j)-h(i)-1)^2; % square trophic coherence of edge i->j
        Isq(j, i) = (h(i)-h(j)-1)^2; % square trophic coherence of edge j->i
        z(i, j) = 1 + exp(gamma * Isq(i, j)); % nomalization factor of edge i->j
        z(j, i) = 1 + exp(gamma * Isq(j, i)); % nomalization factor of edge j->i
        p(i, j) = 1/z(i, j);% probability of edge i->j
        p(j, i) = 1/z(j, i);% probability of edge j->i
        
        %calculate expected number of edges
        n_edge = n_edge + p(i,j) + p(j,i);
    end
end