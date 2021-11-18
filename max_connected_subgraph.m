function Gmax = max_connected_subgraph(G, type)
% max_connected_subgraph  find the largest strongly or weakly connected component

    comp = conncomp(G,'Type',type);
    [~,val] = max(histc(comp,unique(comp)));
    nodes = 1 : numnodes(G);
    Gmax = subgraph(G, nodes(comp==val));
end
