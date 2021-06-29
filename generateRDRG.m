function [A, theta] = generateRDRG(K, m, gamma, a)
% generateRDRG  Generate synthetic network from the directed pRDRG model
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

% initialize model parameters
n = m*K; % total number of nodes
g = 1/K; % parameter g for the Magnetic Laplacian
theta = sort(repmat(linspace(0,2*pi-2*pi/K,K),1,m)+(2*a*rand(1,n)-a));% generate phase angles that form K clusters with additive noice ~unif(-a, a)


% create empty adjacency matrix
A = zeros(n,n);

% preallocate probability matrix
f = zeros(n,n); % probability of edge i<->j
q = zeros(n,n); % probability of edge i->j
l = zeros(n,n); % probability of edge j->i

%loop through all pairs of i and j to update the loglikelihood
for i = 1:(n-1)
    for j = (i+1):n
        
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
        
        
        %update adjacency matrix according to the probability
        rand_number = rand; %draw random number
        if rand_number <= f(i,j)
            A(i,j) = 1; A(j,i) = 1; % generate edge i<->j
        elseif rand_number <= f(i,j)+q(i,j)
            A(i,j) = 1; A(j,i) = 0;  % generate edge i->j
        elseif rand_number <= f(i,j)+q(i,j)+l(i,j)
            A(i,j) = 0; A(j,i) = 1;  % generate edge  j->i
        else
            A(i,j) = 0; A(j,i) = 0; % no edge between i and j
        end

    end
end
end