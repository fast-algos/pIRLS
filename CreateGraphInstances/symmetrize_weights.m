% Function to symmetrize weights (which we always do)
function [knn, wnn, k_neigh] = symmetrize_weights(n, m, k_neigh, knn, wnn)

% Assemble & Symmetrize the weight-matrix:
W = sparse((1:n+m)'*ones(1, k_neigh), knn, wnn); W = 0.5*(W + W');

% Compute Start Index for neighbors of each 'x'.
[row, col] = find(W); % Note: index columns are already sorted, and it's ok
                      % to loop across 'col' because W is symmetric.

% Figure out length of each set of neighbors:
start = [1; find(diff(col)) + 1; length(col)+1]; neigh_x = diff(start);

%% Compute new knn & weights matrices (in place):
knn = ones(n+m, max(neigh_x)); wnn = zeros(n+m, max(neigh_x));
for i = 1 : n+m
    knn(i, 1 : neigh_x(i)) = row(start(i) : start(i+1) - 1);
    wnn(i, 1 : neigh_x(i)) = W(i, knn(i, 1:neigh_x(i)));
end

k_neigh = size(knn, 2); 

end     % End of function to symmetrize weights