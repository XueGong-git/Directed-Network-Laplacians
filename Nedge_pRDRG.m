function [n_edge] = Nedge_pRDRG(theta, gamma, g)
% Nedge_pRDRG   Calculate expected number of edges of the directed pRDRG
%   
% - theta   Estimated phase angles of the nodes
% - gamma   Decay parameter for the directed pRDRG model
% - g       Parameter for the directed pRDRG model

% OUTPUTS
% - n_edge Expected number of edges

n_nodes = length(theta);
n_edge = 0;
f = zeros(n_nodes, n_nodes); % preallocated trophic incoherence
q = zeros(n_nodes, n_nodes); % preallocate normalization factors
l = zeros(n_nodes, n_nodes); % preallocate probabilities of edges
for i = 1 : n_nodes - 1
    for j = i+1 : n_nodes
        
        beta = theta(i) - theta(j); % phase difference
        %calculate terms inside the exponential function
        w1 = 0; %  i<->j
        w2 = gamma*(1-2*cos(beta)+cos(beta+2*pi*g)); %  i->j
        w3 = gamma*(1-2*cos(beta)+cos(beta-2*pi*g)); %  j->i
        w4 = gamma*(2-2*cos(beta)); %  i,j disconnected

        z = exp(w1) + exp(w2)+exp(w3)+exp(w4); % calculate normalization factor
        f(i,j) = 1/z; % calculate the probability of edge i<->j
        q(i,j) = exp(w2)/z; % calculate the probability of edge i->j
        l(i,j) = exp(w3)/z; % calculate the probability of edge j->i
        
        %calculate expected number of edges
        n_edge = n_edge + 2*f(i,j) + q(i,j) + l(i,j);
    end
end