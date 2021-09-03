function [opt_k_ML, Nnodes, Nedges] = optimal_g(input, G_lcc)
% find the optimal parameter g for the Magnetic Laplacian based on the
% maximum likelihood
%
% INPUTS
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
%   G_lcc   Graph object corresponding to the largest connected component
%   of the network
%

% OUTPUTS
% - opt_k_ML  Inverse of the optimal g
% - Nnodes  Number of nodes
% - Nedges Number of edges
%
% DEPENDENCIES
% - meigenmaps -- checked
% - RDRG_likelihood -- checked


%parameters for grid search
kmax = 6; k_test = linspace(2, kmax, kmax-1);
g_ML = 1./k_test; % test likelihood on g_ML
test_gamma = linspace(0,20,41);
Ln = zeros(size(g_ML,2), length(test_gamma));

%load and prepare data
[Nnodes, ~] = size(G_lcc.Nodes);
[Nedges, ~] = size(G_lcc.Edges);
idx_rand = randperm(Nnodes);% shuffle the nodes
G_rand = reordernodes(G_lcc,idx_rand);

for i = 1: size(g_ML,2)
    g = g_ML(i); % parameter g for Magnetic laplacian
    G = meigenmaps(G_rand, g); % solve eigenproblem ofthe Magnetic Laplacian
  
    %calculate log-likelihood
    for j = 1: length(test_gamma)
        Ln(i,j) = RDRG_likelihood(G.Nodes.phase, G, test_gamma(j), g);
    end

end


%find optimal g and gamma 
[idx_g,~] = find(Ln == max(Ln(:)));
%opt_g = g_ML(idx_g);
opt_k_ML = k_test(idx_g);

% plot likelihood of graph over gamma
plot(test_gamma, Ln(1,:), '-*','LineWidth', 1.5);
hold on;
xlabel('\gamma','FontSize', 30)
ylabel('Log-likelihood','FontSize', 30)
plot(test_gamma, Ln(2,:), '--o','LineWidth', 1.5);
plot(test_gamma, Ln(3,:), ':*','LineWidth', 1.5);
plot(test_gamma, Ln(4,:), '-o','LineWidth', 1.5);
plot(test_gamma, Ln(5,:), '-.*','LineWidth', 1.5);
lgd = legend({'g = 1/2','g = 1/3','g = 1/4','g = 1/5','g = 1/6'},'FontSize', 20, 'FontWeight','Bold', 'Location', 'southeast');
set(lgd,'Interpreter','latex');
set(gca,'FontSize',30)
ax = gca;
exportgraphics(ax,strcat('plots/opt_g_input=',num2str(input),'.eps'),'Resolution',300) 
hold off;

end
