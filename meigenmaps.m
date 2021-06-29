function [G] = meigenmaps(G,g)
% meigenmaps Estimate phase angles of nodes using the Magnetic Laplacian
% INPUTS
%   G   graph object of the network
%   g   parameter for the Magnetic Laplacian
%
% OUTPUTS
% - G  Graph object of the network, where G.Nodes.phase is the estimated phase
% angle

    A = adjacency(G,'weighted'); % get adjacency matrix
    Ws = (A + A.')/2;% symmetric weights.
    alpha = sparse(A - A.');  % edge direction.
    d = sum(Ws, 2); Deg = diag(d); % degree matrix
    Tg = exp(1).^(2*pi*1i*g*alpha.');  % Transporter
    Lg = Deg - Ws.*Tg; % magnetic Laplacian.
    [V,~] = eigs(Lg,size(A,1),'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue
    Phi = angle(V); % phases corresponding to all eigenvectors
    G.Nodes.phase = mod(Phi(:,1), 2*pi);% phases corresponding to smallest eigenvalue
    
end
