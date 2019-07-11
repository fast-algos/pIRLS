 % We generate random matrices and run our IRLS algorithm on them. We will
 % set this up to run the algorithm for the unconstrained version.
 % Data dimensions are m x n, where n < m.
 
 n = 800;
 m = 1000;
 eps = 1e-8;
 p = 8;
 A = rand(m,n);
 b = rand(m,1);
 C = zeros(n,n);
 d = zeros(n,n);
 pNorm(eps,A,b,p,C,d);