% Written by Mauricio Flores, last edited 08/24/2018.
% This function computes 'knn' and 'wnn' for a given 'X', 'Y' set.

function [knn, wnn, k_neigh, h] = compute_knn_wnn(X, n, m, k_neigh, local, symm) 

% Apply default inputs:
if nargin < 5; local = 0; end
if nargin < 6; symm = 1; end

% Computing KNN is 99% of the cost:
[knn, dnn] = knnsearch(X, X, 'k', k_neigh);

% Compute length-scale. If 'local' = 1, then we compute a different
% length scale for each vertex in the graph (symmetry is not lost though).
if local == 0
    h = max(dnn(:))/2; wnn = exp(-dnn.^2/h^2); 
else
    h = max(dnn, [], 2)/2; wnn = exp(-(dnn./h).^2);
end

if symm == 1
    [knn, wnn, k_neigh] = symmetrize_weights(n, m, k_neigh, knn, wnn);
end

h = mean(h); % Need to return a single 'h':

end  % End of compute_weights.
