function Gmax = max_connected_subgraph(G, type)
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

    comp = conncomp(G,'Type',type);
    [~,val] = max(histc(comp,unique(comp)));
    nodes = 1 : numnodes(G);
    Gmax = subgraph(G, nodes(comp==val));
end