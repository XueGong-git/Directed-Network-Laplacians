% Generate figures and tables for paper
% DEPENDENCIES
% - load_data
% - optimal_g -- checked
% - model_comparison

clear % clear all variables from workspace
% create file to save results for Table 1
fid = fopen('model_comparison_results','a'); 
fprintf(fid,'%s %s %s %s %s\n', ["Data", "Nnodes", "Nedges", "opt_k", "ln_p_ratio"]); %specify column names

% parameters for synthetic networks
m = 100; % number of nodes per cluster
K = 5; % number of clusters
a = 0.2; % random noise parameter
gamma = 5; % decay parameter
g_RDRG=1/K; % g parameter for pRDRG

for input = 1:15
    [G, ~] = load_data(input,K,m,gamma,a); % get the graph object from the network data
    [opt_k, Nnodes, Nedges] = optimal_g(input,G); % select optimal parameter g for the Magnetic Laplacian
    ln_p_ratio = model_comparison(input,opt_k,G); % compare the liklihood of models
    fprintf(fid,'%d %d %d %d %2d\n', [input, Nnodes, Nedges, opt_k, ln_p_ratio]); %save results to file

end

fclose(fid); %close file
