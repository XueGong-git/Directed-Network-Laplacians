function [G_lcc,A_lcc] = load_data(input, K, m, gamma, a)
% load_data  Load network data, extract the largest connected component, and build the adjacency matrix
% 
% INPUTS
% REQIORED INPUTS
%   input   Integer indicating the input network data
%   input = 1 Synthetic network from directed pRDRG model
%           2 Synthetic network from trophic RDRG model
%           3 Florida Bay food web
%           4 Word adjacency matrix
%           5 Political blogsphere
%           6 Dunnhumby shopping basket data
%           7 Reopen of venue
%           8 US 2015 inflow outflow
%           9 1998 trade network - top 5 export partners
%           10 C-elegans-frontal neural network
%           11 S. cerevisiae transcriptional regulation network
%           12 Transportation reachability
%           13 Influence matrix
%           14 Country to country flight matrix
%           15 US migration
%
% OPTIONAL INPUTS
% - K  Number of clusters, required for synthetic networks when input = 1 or 2
% - n  Number of nodes per cluster, required for synthetic networks when input = 1 or 2
% - gamma  Decay parameter, required for synthetic networks when input = 1 or 2
% - a  Parameter for the additive noise, required when input = 1 or 2
%
% OUTPUTS
% - G_lcc  Graph object corresponding to the largest connected component of the input
% - A_lcc  Adjacency matrix of the largest connected component
% 
% DEPENDENCIES
% - generateRDRG
% - multilevel_model
% - max_connected_subgraph


    if input == 1 %pRDRG model
        [A, theta] = generateRDRG(K, m, gamma, a); %
        G = digraph(A);
        G.Nodes.theta = theta';
        
    elseif input == 2 %trophic RDRG model
        [A, h] = multilevel_model(K, m, gamma, a);
        G = digraph(A);
        G.Nodes.h = h';

    elseif input ==3  %food web, original file can be downloaded from https://snap.stanford.edu/data/Florida-bay.html
        fileID = fopen('datasets/Florida-Bay.txt','r'); 
        formatSpec = '%d %d'; sizeA = [2 Inf]; in = fscanf(fileID,formatSpec,sizeA); edge = in'; % read edges from txt file
        G = digraph(categorical(edge(:,1)),categorical(edge(:,2)));
        NodesList = readtable('datasets/Florida-bay-meta.csv', 'ReadVariableNames',true, 'ReadRowNames',true, 'TextType', 'char');
        NodesList = table2cell(NodesList);
        G.Nodes.Name = NodesList(:,1);
        G = rmnode(G,{'Benthic POC','Water POC',  'DOC'});
        
    elseif input == 4 %word adjacency, original file can be downloaded from http://www-personal.umich.edu/~mejn/netdata/adjnoun.zip
        s = tdfread('datasets/word_adjacency_edges.csv', ';');
        G = digraph(categorical(s.Source),categorical(s.Target));
        
    elseif input == 5 %political blogsphere, original file can be downloaded http://www-personal.umich.edu/~mejn/netdata/polblogs.zip
        s = tdfread('datasets/polblogs_edges.csv', ',');
        G = digraph(categorical(s.Source),categorical(s.Target));
        G = simplify(G);
        
    elseif input == 6 %Dunnhumby 'The Complete Journey' shopping data, original data can be downloaded from https://www.dunnhumby.com/source-files/
        load('datasets/condP');
        threshold = 0.4;
        A = (condP > threshold) - eye(length(condP));
        G = digraph(A);
        G.Nodes.Name = commodity_name;
        
    elseif input == 7 
        %Reopen of venue, the matrix is from Figure 1 in the paper Benzell SG, Collis A, Nicolaides C.Rationing social contact 
        %during the COVID-19 pandemic: Transmission risk and social benefitsof US locations.Proceedings of the National Academy 
        %of Sciences. 2020;117(26):14642–14644.
        %matrix entries with value 1 correspond to the yellow cells in Figure 1
        load('datasets/melted_cormat');
        G = digraph(M, names); 
        
    elseif input == 8 %US 2015 input-output matrix https://stats.oecd.org/Index.aspx?DataSetCode=IOTSI4_2018
        T = readtable('datasets/US_io_2015.csv', 'ReadRowNames', 1);
        W = table2array(T); % input weight matrix
        A = zeros(size(W));
        thereshold = prctile(W,70, 'all'); A(W >= thereshold) = 1; A(W == W') = 0;%binarize the matrix 
        G = digraph(A);
        G.Nodes.Name = T.Properties.RowNames;  
        
    elseif input == 9 %world trade relations http://www.economicswebinstitute.org/worldtrade.htm
        T = readtable('datasets/trade_1998.csv', 'ReadRowNames', 1);
        W = table2array(T);
        G = digraph(W);
        G.Nodes.Name = T.Properties.RowNames;
        
    elseif input == 10 %C-elegans-frontal, original data downloaded from http://snap.stanford.edu/data/C-elegans-frontal.html
        fileID = fopen('datasets/C-elegans-frontal.txt','r');
        formatSpec = '%d %d';
        sizeA = [2 Inf];
        in = fscanf(fileID,formatSpec,sizeA);
        edge = in';
        G = digraph(categorical(edge(:,1)),categorical(edge(:,2)));
        
    elseif input == 11 %S-cerevisiae, original data downloaded from http://snap.stanford.edu/data/S-cerevisiae.html
        fileID = fopen('datasets/S-cerevisiae.txt','r');
        formatSpec = '%d %d %d';
        sizeA = [3 Inf];
        in = fscanf(fileID,formatSpec,sizeA);
        edge = in';
        G = digraph(categorical(edge(:,1)),categorical(edge(:,2)));
        
    elseif input == 12 %transportation reachability, original data downloaded from http://snap.stanford.edu/data/reachability.html
        fileID = fopen('datasets/reachability.txt','r');
        formatSpec = '%d %d %d';
        sizeA = [3 Inf];
        in = fscanf(fileID,formatSpec,sizeA);
        edge = in';
        G = digraph(categorical(edge(:,1)),categorical(edge(:,2)));
        
    elseif input == 13 %influence matrix, from Cole A.The influence matrix methodology: A technical report.Landcare Research Contract Report: LC0506/175. 2006.
        T = readtable('datasets/imatrix.csv', 'ReadRowNames', 1);
        W = table2array(T);
        W = W - diag(diag(W));
        A = select_top_partners(W, 1, 2); 
        G = digraph(A);
        G.Nodes.Name = T.Properties.RowNames;  
        imagesc(A);
        
    elseif input == 14 %country to country flight matirx, original data downloaded from https://www.visualizing.org/global-flights-network/
        T = readtable('datasets/airport_CnToCn_ajc.csv', 'ReadRowNames', 1);
        W = table2array(T);
        A = select_top_partners(W,10, 2);
        G = digraph(A);
        G.Nodes.Name = T.Properties.RowNames; 
    
    elseif input == 15 %US migration, original data downloaded from https://www.census.gov/content/census/en/library/publications/2003/dec/censr-8.html
        W = readtable('datasets/US state migration.csv', 'ReadVariableNames',true, 'ReadRowNames',true, 'Format','auto');
        NodesList = readtable('datasets/US_states_list.csv', 'ReadVariableNames',false, 'ReadRowNames',false, 'TextType', 'char');
        W = W{:,:}; W=W'; 
        A = zeros(size(W));
        threshold = prctile(W,70, 'all'); A( W >= threshold) = 1; A(W == W') = 0; %binarize the matrix 
        G = digraph(A);
        NodesList = table2cell(NodesList);
        G.Nodes.label = NodesList;

    end
    
    if any([2, 4, 6, 7, 11]== input)
            G_lcc = max_connected_subgraph(G,'weak'); % extract largest weakly connected component
    else
         G_lcc = max_connected_subgraph(G,'strong'); % extract largest strongly connected component
    end
    A_lcc = sparse(adjacency(G_lcc, 'weighted')); % output adjacency matrix as sparse matrix

end

