function [M] = updateM(X,Xbar,M,F,N,k,p,r)
%
% serves to use multiplicative update rules to update nonnegative tensor M (equation 4.17 in paper Mixtures of Block Models for Brain Networks, SDM 2018)
% (By Zilong Bai, KDD Lab @ University of California, Davis)
% Input
%	X: N x N x r tensor of multigraph. Each frontal slice X(:,:,j) is the affinity/similarity matrix of one graph.
%       Xbar: N x N x r tensor of multigraph reconstructed with current F and M
%       F: N x k x p tensor. N is the number of vertices from the original graphs, k is the number of blocks in each F(:,:,i), p is the number of latent block models.
%       M: k x k x p x r tensor. r is the number of graphs in the multigraph X.
%       N: the number of vertices of each single graph
%       k: the number of blocks in each block structure
%       p: the number of different block structures to find
%       r: the number of single graphs in the multigraph X
% Output
%	M: k x k x p x r tensor. Updated M.

epsilon = 1e-6; % ignorable positive value to avoid issues in multiplicative update rules caused by intermediate zero elements.

for i = 1:p
for j = 1:r

numerator = F(:,:,i)'*X(:,:,j)*F(:,:,i);

denominator = F(:,:,i)'*Xbar(:,:,j)*F(:,:,i);

denominator = denominator + ones(size(denominator)).*epsilon;

nd = numerator./denominator;

nd(find(nd<=0)) = epsilon;

nd(find(isnan(nd)==1)) = epsilon;

tempMij = M(:,:,i,j).*nd;

tempMij(find(isnan(tempMij)==1)) = epsilon;

M(:,:,i,j) = tempMij;

end
end
