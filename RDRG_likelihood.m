function [Ln] = RDRG_likelihood(theta, G, gamma, g)
% RDRG_likelihood  Calculate the likelihood that the network came from the directed pRDRG model
% INPUTS
%   
% - theta   Estimated phase angles of the nodes
% - G       Graph object of the network
% - gamma   Decay parameter for the directed pRDRG model
% - g       Parameter for the directed pRDRG model
%
% OUTPUTS
% - Ln lliklihood of directed pRDRG model
% 


[Nnodes, ~] = size(theta); %number of nodes
Ln = 0;
A = adjacency(G);

for i = 1:Nnodes-1
    for j = i+1:Nnodes
        
        beta = theta(i) - theta(j);% phase difference
        %calculate terms inside the exponential function
        w1 = 0;%  i<->j
        w2 = gamma*(1-2*cos(beta)+cos(beta+2*pi*g));%  i->j
        w3 = gamma*(1-2*cos(beta)+cos(beta-2*pi*g));%  j->i
        w4 = gamma*(2-2*cos(beta));%  i,j disconnected
        lnZ = log(exp(w1) + exp(w2) + exp(w3) + exp(w4)); %normalization factor
        
        %update log-likelihood of the graph by adding the log-liklihood of
        %edges
        if A(i,j)==1 && A(j,i)==1%  i<->j
            Ln = Ln - lnZ;
            
        elseif A(i,j)==1 && A(j,i)==0%  i->j
            Ln = Ln + w2 - lnZ;
         
        elseif A(i,j)==0 && A(j,i)==1%  j->i
            Ln = Ln + w3 - lnZ;
            
        else%  i,j disconnected
            Ln = Ln + w4 - lnZ;
            
        end
        
    end
               
end