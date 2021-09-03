function [ln_p_ratio] = model_comparison(input, opt_k, G)
% model_comparison  Calculate and compare the likelihood that the given network is generated from the directed pRDRG model, and the trophic RDRG model
%
% INPUTS
%   - input   Integer indicating the input network data
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
%   - opt_k   The inverse of the optimal g for the Magnetic Laplacian
%   - G   Graph object corresponding to the largest connected component
%         of the network
%

% OUTPUTS
% - ln_p_ratio log of the ratio between the liklihood of directed pRDRG and the
% trohpic RDRG model
%
% DEPENDENCIES
% - meigenmaps -- checked
% - RDRG_likelihood -- checked
% - multilevel_likelihood -- checked
% - Nedge_pRDRG -- checked
% - Nedge_TL -- checked
% - levels (TC Toolbox)

g = 1/opt_k; % get the optimal g for the Magnetic Laplacian
[Nnodes, ~] = size(G.Nodes);
A = adjacency(G);
idx_rand = randperm(Nnodes);% shuffle the nodes
G_rand = reordernodes(G,idx_rand);
A_rand = A(idx_rand,idx_rand);

%inverse the node attributes
G_ML = meigenmaps(G_rand, g); % estimage phase angles using the Magnetic Laplacian
[h_est] = levels(A_rand); % estimate trophic level using the Trophic Laplacian

%preallocate arrays
test_gamma = linspace(0,20,41);% test points for gamma
Ln_ML = zeros(1, length(test_gamma));% log-likelihood of directed pRDRG model
Ln_TL = zeros(1, length(test_gamma));% log-likelihood of Trohpic RDRG model
nedge_exp_ML = zeros(1, length(test_gamma));% expected number of edges of directed pRDRG 
nedge_exp_TL = zeros(1, length(test_gamma));% expected number of edges of trophic RDRG 

for k = 1: length(test_gamma)
    gamma = test_gamma(k);
    Ln_ML(k) = RDRG_likelihood(G_ML.Nodes.phase, G_rand, gamma, g); %likelihood of pRDRG
    Ln_TL(k) = multilevel_likelihood(h_est, G_rand, gamma); %likelihood of trophic level
    nedge_exp_ML(k) = Nedge_pRDRG(G_ML.Nodes.phase, gamma, g); %expected number of edges of directed pRDRG
    nedge_exp_TL(k) = Nedge_TL(h_est, gamma); %expected number of edges of trophic pRDRG model

end

%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot original matrix and the attributes
img = imagesc(A,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2))); 
set(gca,'FontSize',30)
%img.LineWidth = 2;
ax = gca;
exportgraphics(ax,strcat('plots/original_input',num2str(input),'.eps'),'Resolution',300) 

% plot phase angles estimated using the Magnetic Laplacian; and the trophic
% levels estimated using the Trohpic Laplacian
if Nnodes <= 100 %plot network is the number of nodes is smaller than 100
    
    % plot phase angles estimated using the Magnetic Laplacian
    if input == 3 % Florida bay food web
        fig = plot(G_ML,'XData', cos(G_ML.Nodes.phase) ,'YData', sin(G_ML.Nodes.phase), 'NodeLabel', G_ML.Nodes.Name, 'ArrowPosition', 0.8, 'ArrowSize', 10); %plot network on polar coordinate
        text(fig.XData+[0 0.04 0 0 0 -0.23 0.3 0 0 0 0 -0.6 ], ...
             fig.YData+[-0.12 0.04 -0.05 0 0 -0.2 -0.05 0 0 -0.17 0 -0.15 ],fig.NodeLabel, ...
            'VerticalAlignment','Bottom',...
             'HorizontalAlignment', 'left',...
              'FontSize', 20)% adjust node label position 
        fig.NodeLabel = {}; % Remove old labels
        highlight(fig,'Edges',[1, 12, 24, 25],'EdgeColor', 'r', 'LineStyle' , '--') %highlight edges in the clockwise direction in red
    
    elseif input == 14 % influence matrix
        fig = plot(G_ML,'XData', cos(G_ML.Nodes.phase) ,'YData', sin(G_ML.Nodes.phase), 'NodeLabel', G_ML.Nodes.Name, 'ArrowPosition', 0.8, 'ArrowSize', 10); %plot on polar coordinate
        text(fig.XData+[0 -1.2 -0.15 -0.1 0 0 -0.25 0 0 0.1 0 0 -0.6 0 ], ...
             fig.YData+[0.11 0.05 0 -0.15 0 0 -0.25 0 0 -0.05 0.05 0.02 0.05 -0.15 ],fig.NodeLabel, ...
            'VerticalAlignment','Bottom',...
             'HorizontalAlignment', 'left',...
              'FontSize', 15);% adjust node label position 
        fig.NodeLabel = {}; % remove old labels
        highlight(fig,'Edges',[4, 14, 25, 26, 29, 32], 'EdgeColor','r' ); %highlight downward edges in red
        highlight(fig, [1, 2, 7, 8, 9, 10], 'NodeColor', 'r'	) %highlight ecological factors in red
        highlight(fig, [5], 'NodeColor',[0.4940 0.1840 0.5560]	) %highlight policy in purple      
   
    else     
        fig = plot(G_ML,'XData', cos(G_ML.Nodes.phase) ,'YData', sin(G_ML.Nodes.phase), 'ArrowPosition', 0.8); %plot on polar coordinate
    end
    
    % customize plot
    xlim([-1.2, 1.8]); ylim([-1.2, 1.2]);
    fig.MarkerSize = 7; fig.NodeFontSize = 15; fig.LineWidth = 2; fig.LineStyle = '-';
    ax = gca; ax.XAxisLocation = 'origin';ax.YAxisLocation = 'origin';
    xticks([-1  1]); xticklabels({'\pi','0'});
    yticks([-1  1]); yticklabels({'3\pi/2','\pi/2'});
    ax.XRuler.TickLabelGapOffset = 10; ax.YRuler.TickLabelGapOffset = 35;    % move the ticklabel
    label_y = ylabel('Phase angle','FontSize', 13);
    label_y.Position(1) = label_y.Position(1) - 1.7; % change horizontal position of ylabel
    label_y.Position(2) = label_y.Position(2) - 0.4; % change vertical position of ylabel
    set(label_y, 'Rotation', 90,'VerticalAlignment','middle', 'HorizontalAlignment','right')
    set(gca,'FontSize',30) ;
    exportgraphics(ax,strcat('plots/ML_eigmap_input=', num2str(input) ,'_g=',num2str(round(g,2)),'.eps'),'Resolution',300) 

    % plot trophic levels estimated using the Trophic Laplacian
    p=plot(G_rand,'layout','layered');
    xdata=get(p,'XData'); 
    close
    
    if input == 3% Florida bay food web
        fig=plot(G_rand,'XData',xdata,'YData',h_est ,'NodeLabel', G_rand.Nodes.Name, 'ArrowPosition', 0.8, 'ArrowSize', 10);
        text(fig.XData+[0 0 -1 -1 0 -0.5 -3.5 -0.8 0 -0.3 0 -0.5 ], ...
             fig.YData+[0 0 0 -0.1 0 0.03 -0.1 0.05 -0.13 0 0 0 ],fig.NodeLabel, ...
            'VerticalAlignment','Bottom',...
             'HorizontalAlignment', 'left',...
              'FontSize', 20)% adjust node label position
        highlight(fig,'Edges',[5, 11, 13, 21, 24, 28], 'EdgeColor','r' );
        fig.NodeLabel = {};   % Remove old labels
        
    elseif input == 14% influence matrix
        fig=plot(G_rand,'XData',xdata,'YData',h_est ,'NodeLabel', G_rand.Nodes.Name, 'ArrowPosition', 0.8, 'ArrowSize', 10);
        text(fig.XData+[-0.5 -2 -1.5 -1.2 -3.2 0 0 -2 0 -2 0 0 -2.5 0 ], ...
             fig.YData+[-0.1 0.05 0.05 -0.1 0.05 0 0.05 0 0.05 -0.2 0 0 0 0 ],fig.NodeLabel, ...
            'VerticalAlignment','Bottom',...
             'HorizontalAlignment', 'left',...
              'FontSize', 15)% adjust node label position
        highlight(fig,'Edges',[5, 6, 17, 19, 21, 23, 25, 32], 'EdgeColor','r' );
        highlight(fig, [1, 2, 7, 8, 9, 10], 'NodeColor', 'r'	) %highlight ecological factors in red
        highlight(fig, [5], 'NodeColor',[0.4940 0.1840 0.5560]	)  %highlight policy
        fig.NodeLabel = {};   % Remove old labels
    else
        fig=plot(G_rand,'XData',xdata,'YData',h_est, 'ArrowPosition', 0.8);
    end

    % customize plot
    fig.MarkerSize = 6;
    fig.NodeFontSize = 12;
    fig.LineWidth = 2;
    ylabel('Trophic level','FontSize', 13);
    set(gca,'xticklabel',{[]})
    fig.LineStyle = '-';
    set(gca,'fontsize',30);
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/TL_eigmap_input=', num2str(input),'_g=',num2str(round(g,2)),'.eps'),'Resolution',300) 

end


% plot the matrix reordered using the phase angle 
[~,idx_ML] = sort(mod(G_ML.Nodes.phase, 2*pi));
W_ML = A_rand(idx_ML,idx_ML);
imagesc(W_ML,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) 
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/reordered_ML_input_', num2str(input) ,'_g=',num2str(round(g,2)),'.eps'),'Resolution',300) 


% plot the matrix reordered using the Trophic Level
[~,idx_TL] = sort(h_est);
A_TL = A_rand(idx_TL,idx_TL);
imagesc(A_TL,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/reordered_TL_input_', num2str(input) ,'_g=',num2str(round(g,2)),'.eps'),'Resolution',300) 

% plot estimated phase angle vs. truth for the syntheitc network from the
% directed pRDRG model
if input == 1
    plt = scatter(G_ML.Nodes.theta, G_ML.Nodes.phase, '');
    xlabel('True \theta','FontSize', 13);
    ylabel('Estimated \theta','FontSize', 13);
    set(gca,'fontsize',30);
    plt.LineWidth = 2;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/theta_est_input_', num2str(input) ,'_g=',num2str(round(g,2)),'.eps'),'Resolution',300) 
end

% plot estimated trophic level vs. truth for the syntheitc network from the
% Trophic RDRG
if input == 2
    plt = scatter(G_rand.Nodes.h, h_est, '');
    xlabel('True trophic level','FontSize', 13);
    ylabel('Estimated trophic level','FontSize', 13);
    set(gca,'fontsize',30);
    plt.LineWidth = 2;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/h_est_input_', num2str(input) ,'_g=',num2str(round(g,2)),'.eps'),'Resolution',300) 
end



% highlight gamma estimators
%point estimate of gamma
[~, mean_edge_idx] = min(abs(nedge_exp_ML - numedges(G_ML)));
gamma_mean_edge = test_gamma(mean_edge_idx); % point estimate of gamma for the directed pRDRG model
[~, mean_edge_idx_TL] = min(abs(nedge_exp_TL - numedges(G_ML)));
gamma_mean_edge_TL = test_gamma(mean_edge_idx_TL); % point estimate of gamma for the Trophic RDRG model

%maximum likelihood of gamma
[~, max_Ln_ML_idx] = max(Ln_ML);
gamma_max_Ln_ML = test_gamma(max_Ln_ML_idx); 
[~, max_Ln_TL_idx] = max(Ln_TL);
gamma_max_Ln_TL = test_gamma(max_Ln_TL_idx); 

% plot likelihood
plt = plot(test_gamma, Ln_ML, '', 'LineWidth',1.5);
hold on;
plot(test_gamma, Ln_TL, '--r', 'LineWidth',1.5);
plot(gamma_mean_edge, Ln_ML(mean_edge_idx), '*k', 'MarkerSize',10, 'LineWidth',2);
plot(gamma_max_Ln_ML, Ln_ML(max_Ln_ML_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
plot(gamma_mean_edge_TL, Ln_TL(mean_edge_idx_TL), '*r', 'MarkerSize',10, 'LineWidth',2);
plot(gamma_max_Ln_TL, Ln_TL(max_Ln_TL_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
legend({'Magnetic Laplacian','Trophic Level', 'Point estimate', 'MLE'},'FontSize', 20,'Location','southwest');
xlabel('\gamma','FontSize', 13);
ylabel('Log-likelihood','FontSize', 13);
set(gca,'fontsize',30);
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,strcat('plots/model_comparison_input_', num2str(input), '_g=', num2str(round(g,2)),'.eps'),'Resolution',300) 
hold off;


% calculate difference between the log likelihood (or log of the likelihood
% ratio)
ln_p_ratio = Ln_ML(max_Ln_ML_idx) - Ln_TL(max_Ln_TL_idx);

end
