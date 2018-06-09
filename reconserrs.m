function [err] = reconserrs(X,Xbar)
%
% serves to calculate the relative reconstruction error of the multigraph.
% 
% Input
%	X: N x N x r tensor of multigraph. Each frontal slice X(:,:,j) is the affinity/similarity matrix of one graph.
%	Xbar: N x N x r tensor of multigraph reconstructed with F and M before the current iteration of updating.
% Output
%	err: the mean relative reconstruction error of all the single graphs that comprise the multigraph X.

[N,N1,r] = size(X);

diffnorm = zeros(r,1);
orinorm = zeros(r,1);

for j = 1:r
diffnorm(j) = norm((X(:,:,j)-Xbar(:,:,j)),'fro');
orinorm(j) = norm(X(:,:,j),'fro'); 
end

relerr = diffnorm./orinorm;

err = mean(relerr);
