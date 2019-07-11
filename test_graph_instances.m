 % We use the graph instances generate by Mauricio Flores et al. and run
 % our IRLS algorithm to find the minimum $\ell_p$ laplacian.
 
 eps = 1e-8;
 n = 1000; 
 p = 8;
 nearest_neighbours = 10;
 [A,b,C,d] = create_graph_matrices(n,p,nearest_neighbours);
 pNorm(eps,A,b,p,C,d);