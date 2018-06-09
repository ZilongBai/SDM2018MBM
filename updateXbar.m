function [Xbar,Xbarij] = updateXbar(X,M,F)
%
% serves to calculate intermediate variables, in particular graphs reconstructed with current F and M, to reduce computationl redundancy in implementing update rules for F and M. (By Zilong Bai, KDD Lab @ University of California, Davis)
% Input
%	X: N x N x r tensor of multigraph. Each frontal slice X(:,:,j) is the affinity/similarity matrix of one graph.
%	M: k x k x p x r current tensor of sets of weighted mixing matrices. 
% 	F: N x k x p current tensor of block structures. 
%
% Output
%	Xbar: N x N x r tensor of multigraph reconstructed with current F and M
%	Xbarij: N x N x p x r tensor of intermediate reconstructed graphs. Xbarij(i,j) reconstructed with F(:,:,i) and M(:,:,i,j).
 
[N,N1,r] = size(X);
[N2,k,p] = size(F);

Xbar = zeros(N,N,r);

Xbarij = zeros(N,N,p,r);

for j = 1:r

temp = zeros(N,N,p);

for i = 1:p

temp(:,:,i) = F(:,:,i)*M(:,:,i,j)*F(:,:,i)';

Xbarij(:,:,i,j) = temp(:,:,i);

end

Xbar(:,:,j) = tensorvectordot(temp,ones(p,1));

end
