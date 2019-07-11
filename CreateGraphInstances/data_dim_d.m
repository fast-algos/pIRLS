function [X, Y, g, m] = data_dim_d(n, d)

% Written by Mauricio Flores, last edited on 08/22/18.

% We have decided to use three types of datasets, as follows
% (1) Two dimensional, uniform, with two labels.
% (2) High-dimensional, uniform, with few labels.
% (3) Two clusters, with different functions inside.

% This corresponds to dataset (2) within that list.
% We generate a uniform distribution in [0,1]^d, 
% with labeled points & labels randomly chosen.

% Unlabeled points
X = rand(n, d); 

% Labeled points
m = 10;
g = rand(m, 1);
Y = rand(m, d);

end