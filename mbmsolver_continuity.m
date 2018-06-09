function [F,M,errs,Xbar] = mbmsolver_continuity(X, k, p, MAX, pc, theta)
%
% solves the Block Structure Set Linked Matrix tri-Factorization with Spatial Continuity Regularization (BSSLMtF-SCR) for sets of  weighted mixing matrices M and small set of block structures F (for Algorithm 1 Mixtures of Block Models Discovery in Paper Mixtures of Block Models for Brain Networks, SDM 2018) 
% (by Zilong Bai, KDD Lab @ University of California, Davis)
%
% The small set of block structures F is solved with multiplicative update rules in updateF_continuity.m
% The sets of weighted mixing matrices M are solved with multiplicative update rules in updateM.m 
%% The detailed formulation and deductions can be found in our paper Mixtures of Block Models for Brain Networks, SIAM Data Mining 2018 (SDM'18). Please refer to the related work for techniques involved in developing our method.
%
%% Input
%	X: N x N x r tensor of multigraph. Each frontal slice X(:,:,j) is the affinity/similarity matrix of one graph.
%	k: The number of blocks in each block structure
%	p: The number of different block structures to find
%	MAX: maximum number of iterations
%	pc: nonnegative scalar, Weight of the continuity regularization term, i.e., beta in the paper.
%	Theta: N x N spatial closeness penalty matrix. In particular for the applciation of our method to networks constructed with brain imaging data in this paper, the penalty matrix Theta is calculated according to the spatial coordinates of the centroids of the anatomical regions in the brain.
%
% Output
%	F: N x k x p, the tensor of block structures. F(:,:,i) is the block structure for the i-th block model in mdoeling each graph.
%	M: k x k x p x r, the tensor of sets of weighted mixing matrices. M(:,:,i,j) records the interaction between blocks of the block structure F(:,:,i) for modeling graph X(:,:,j). 
%	errs: vector for monitoring how well the discovered block models fit the data.

[N, N1, r] = size(X);

% Initialization

[F,M,alpha] = InitFMAlpha(N,k,p,r);

for i = 1:MAX

% Calculating intermediate variables to reduce redundancy
[Xbar,Xbarij] = updateXbar(X,M,F);

% Updating F with multiplicative update rule

[F] = updateF_continuity(X,Xbar,M,F, N,k,p,r, pc, theta);

% Update M with multiplicative update rule

[M] = updateM(X,Xbar,M,F,N,k,p,r);

% Record relative reconstruction errors.

errs(i) = reconserrs(X,Xbar);

end

